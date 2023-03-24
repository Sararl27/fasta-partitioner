import os
import pathlib
import re
import lithops


class FastaPartitioner:

    def __init__(self, storage, bucket, key, workers):
        self.storage = storage
        self.bucket = bucket

        self.__generate_fasta_index(key, workers)

    def __get_length(self, min_range, content, data, start_base, end_base, id):
        start_base -= min_range
        end_base -= min_range
        len_base = len(data[start_base:end_base].replace('\n', ''))
        # name_id offset_head offset_bases ->
        # name_id offset_head offset_bases len_bases id
        content[-1] = f'{content[-1]} {len_base} {str(id)}'

    # Generate metadata from fasta file
    def __generate_chunks(self, id, key, chunk_size, obj_size, partitions):
        min_range = id * chunk_size
        max_range = int(fasta_size) if id == num_chunks - 1 else (id + 1) * chunk_size
        data = storage.get_object(bucket=fasta_path.bucket, key=fasta_path.key,
                                  extra_get_args={'Range': f'bytes={min_range}-{max_range - 1}'}).decode('utf-8')
        content = []
        # If it were '>' it would also find the ones inside the head information
        ini_heads = list(re.finditer(r"\n>", data))
        heads = list(re.finditer(r">.+\n", data))

        if ini_heads or data[0] == '>':  # If the list is not empty or there is > in the first byte
            first_sequence = True
            prev = -1
            for m in heads:
                start = min_range + m.start()
                end = min_range + m.end()
                if first_sequence:
                    first_sequence = False
                    if id > 0 and start - 1 > min_range:
                        # If it is not the worker of the first part of the file and in addition it
                        # turns out that the partition begins in the middle of the base of a sequence.
                        # (start-1): avoid having a split sequence in the index that only has '\n'.
                        match_text = list(re.finditer('.*\n', data[0:m.start()]))
                        if match_text and len(match_text) > 1:
                            text = match_text[0].group().split(' ')[0].replace('\n', '')
                            offset = match_text[1].start() + min_range
                            # >> offset_head offset_bases_split ^first_line_before_space_or_\n^
                            content.append(f">> <Y> {str(offset)} ^{text}^")  # Split sequences
                        else:
                            # When the first header found is false, when in a split stream there is a split header
                            # that has a '>' inside (ex: >tr|...o-alpha-(1->5)-L-e...\n)
                            first_sequence = True
                if prev != start:  # When if the current sequence base is not empty
                    # name_id offset_head offset_bases
                    id_name = m.group().replace('\n', '').split(' ')[0].replace('>', '')
                    content.append(f"{id_name} {str(start)} {str(end)}")
                prev = end

            # Check if the last head of the current one is cut. (ini_heads[-1].start() + 1): ignore '\n'
            if len(heads) != 0 and len(ini_heads) != 0 and ini_heads[-1].start() + 1 > heads[-1].start():
                last_seq_start = ini_heads[-1].start() + min_range + 1  # (... + 1): ignore '\n'
                text = data[last_seq_start - min_range::]
                # [<->|<_>]name_id_split offset_head
                # if '<->' there is all id
                content.append(f"{'<-' if ' ' in text else '<_'}{text.split(' ')[0]} {str(last_seq_start)}")
        return content

    def __reduce_generate_chunks(self, results):
         if len(results) > 1:
            results = list(filter(None, results))
            for i, list_seq in enumerate(results):
                if i > 0:
                    list_prev = results[i - 1]
                    # If it is not empty the current and previous dictionary
                    if list_prev and list_seq:
                        param = list_seq[0].split(' ')
                        seq_prev = list_prev[-1]
                        param_seq_prev = seq_prev.split(' ')

                        # If the first sequence is split
                        if '>>' in list_seq[0]:
                            if '<->' in seq_prev or '<_>' in seq_prev:
                                # If the split was after a space, then there is all id
                                if '<->' in seq_prev:
                                    name_id = param_seq_prev[0].replace('<->', '')
                                else:
                                    name_id = param_seq_prev[0].replace('<_>', '') + param[3].replace('^', '')
                                list_seq[0] = rename_sequence(list_seq[0], param, name_id, param_seq_prev[1], param[2])
                            else:
                                list_seq[0] = seq_prev
                            # Remove previous sequence
                            list_prev.pop()
        return list(filter(None, results))

    def __generate_index_file(self, data, file_name):
        output_path = './output_data/'
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        with open(f'{output_path}{file_name}_index.fai', 'w') as f:
            for list_seq in data:
                for sequence in list_seq:
                    f.write(f'{sequence}\n')


    def __generate_fasta_index(self, key, workers):
        fexec = lithops.FunctionExecutor(max_workers=2000, runtime_memory=4096)  # log_level='DEBUG

        # location = obj.split('//')[1].split('/')  # obj.split('/')
        # for_head = location[1:]
        # for_head = '/'.join(for_head)
        # data_bucket_name = location[0]
        fasta = self.storage.head_object(self.bucket, key)
        chunk_size = int(int(fasta['content-length']) / workers)
        # seq_name = location[-1].split('.')[0]

        print('===================================================================================')
        print('metadata chunks: ' + str(chunk_size))
        print('bucket to access data: ' + str(self.bucket))
        print('reference file name: ' + pathlib.Path(key).stem)
        print('fasta size: ' + str(fasta['content-length']) + ' bytes')
        print('===================================================================================')

        map_iterdata = [{'key': key} for _ in range(workers)]
        extra_args = {'chunk_size': chunk_size, 'obj_size': fasta['content-length'], 'partitions': workers}
        fexec.map_reduce(map_function=self.__generate_chunks, map_iterdata=map_iterdata, extra_args=extra_args,
                         reduce_function=self.__reduce_generate_chunks)
        results = fexec.get_result()

        self.__generate_index_file(results, f'{pathlib.Path(key).stem}')

        fexec.clean()

        print('... Done, generated index')
