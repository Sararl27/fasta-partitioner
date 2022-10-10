import os
import pathlib
import re
import lithops


class FastaPartitioner:

    def __init__(self, storage, bucket, key, workers):
        self.storage = storage
        self.bucket = bucket

        self.__generate_fasta_index(key, workers)

    def __get_length(self, min_range, content, data, start_base, end_base):
        start_base -= min_range
        end_base -= min_range
        len_base = len(data[start_base:end_base].replace('\n', ''))
        # name_id num_chunks_has_divided offset_head offset_bases ->
        # name_id num_chunks_has_divided offset_head offset_bases len_bases
        content[-1] = f'{content[-1]} {len_base}'

    # Generate metadata from fasta file
    def __generate_chunks(self, id, key, chunk_size, obj_size, partitions):
        min_range = id * chunk_size
        max_range = int(obj_size) if id == partitions - 1 else (id + 1) * chunk_size
        data = self.storage.get_object(bucket=self.bucket, key=key,
                                       extra_get_args={'Range': f'bytes={min_range}-{max_range - 1}'}).decode('utf-8')
        content = []
        ini_heads = list(
            re.finditer(r"\n>", data))  # If it were '>' it would also find the ones inside the head information
        heads = list(re.finditer(r">.+\n", data))

        if ini_heads or data[
            0] == '>':  # If the list is not empty or there is > in the first byte (if is empty it will return an empty list)
            first_sequence = True
            prev = -1
            for m in heads:
                start = min_range + m.start()
                end = min_range + m.end()
                if first_sequence:
                    first_sequence = False
                    if id > 0 and start - 1 > min_range:  # If it is not the worker of the first part of the file and in addition it
                        # turns out that the partition begins in the middle of the base of a sequence.
                        # (start-1): avoid having a split sequence in the index that only has '\n'
                        match_text = list(re.finditer('.*\n', data[0:m.start()]))
                        if match_text:
                            text = match_text[0].group().split(' ')[0]
                            length_0 = len(data[match_text[0].start():m.start()].replace('\n', ''))
                            offset_0 = match_text[0].start() + min_range
                            if len(match_text) > 1:
                                offset_1 = match_text[1].start() + min_range
                                length_1 = len(data[match_text[1].start():m.start()].replace('\n', ''))
                                length_base = f"{length_0}-{length_1}"
                                offset = f"{offset_0}-{offset_1}"
                            else:
                                length_base = f"{length_0}"
                                offset = f'{offset_0}'
                            # >> offset_head offset_bases_split length/s first_line_before_space_or_\n
                            content.append(f">> <Y> {str(offset)} {length_base} ^{text}^")  # Split sequences
                        else:  # When the first header found is false, when in a split stream there is a split header that has a '>' inside (ex: >tr|...o-alpha-(1->5)-L-e...\n)
                            first_sequence = True
                            start = end = -1  # Avoid entering the following condition
                if prev != start:  # When if the current sequence base is not empty
                    if prev != -1:
                        self.__get_length(min_range, content, data, prev, start)
                    # name_id offset_head offset_bases
                    id = m.group().replace('\n', '').split(' ')[0].replace('>', '')
                    content.append(f"{id} {str(start)} {str(end)}")
                prev = end
            
            if len(heads) != 0 and len(ini_heads) != 0:
                if ini_heads[-1].start() + 1 > heads[
                    -1].start():  # Check if the last head of the current one is cut. (ini_heads[-1].start() + 1): ignore '\n'
                    last_seq_start = ini_heads[-1].start() + min_range + 1  # (... + 1): ignore '\n'
                    self.__get_length(min_range, content, data, prev, last_seq_start) # Add length of bases to last sequence
                    text = data[last_seq_start - min_range::]
                    # [<->|<_>]name_id_split offset_head
                    content.append(
                        f"{'<-' if ' ' in text else '<_'}{text.split(' ')[0]} {str(last_seq_start)}")  # if '<->' there is all id
                else:  # Add length of bases to last sequence
                    self.__get_length(min_range, content, data, prev, max_range)

        return content

    def __reduce_generate_chunks(self, results):
        if len(results) > 1:
            results = list(filter(None, results))
            for i, list_seq in enumerate(results):
                if i > 0:
                    list_prev = results[i - 1]
                    if list_prev and list_seq and '>>' in list_seq[
                        0]:  # If i > 0 and not empty the current and previous dictionary and the first sequence is split
                        param = list_seq[0].split(' ')
                        seq_prev = list_prev[-1]
                        param_seq_prev = seq_prev.split(' ')
                        if '<->' in seq_prev or '<_>' in seq_prev:
                            if '<->' in list_prev[-1]:  # If the split was after a space, then there is all id
                                name_id = param_seq_prev[0].replace('<->', '')
                            else:
                                name_id = param_seq_prev[0].replace('<_>', '') + param[4].replace('^', '')
                            length = param[3].split('-')[1]
                            offset_head = param_seq_prev[1]
                            offset_base = param[2].split('-')[1]
                            list_prev.pop()  # Remove previous sequence
                        else:
                            length = param[3].split('-')[0]
                            name_id = param_seq_prev[0]
                            offset_head = param_seq_prev[1]
                            offset_base = param[2].split('-')[0]
                        list_seq[0] = list_seq[0].replace(f' {param[4]}', '')  # Remove 4rt param
                        list_seq[0] = list_seq[0].replace(f' {param[2]} ',
                                                      f' {offset_base} ')  # [offset_base_0-offset_base_1|offset_base] -> offset_base
                        list_seq[0] = list_seq[0].replace(f' {param[3]}', f' {length}')  # [length_0-length_1|length] -> length
                        list_seq[0] = list_seq[0].replace(' <Y> ', f' {offset_head} ')  # Y --> offset_head
                        list_seq[0] = list_seq[0].replace('>> ', f'{name_id} ')  # '>>' -> name_id
        return results

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

        # fexec.plot()
        fexec.clean()

        print('... Done, generated index')


class FunctionsFastaIndex:
    def __init__(self, path_index_file):
        self.data_path = path_index_file

    def get_info_sequence(self, identifier):
        length = offset_head = offset = None
        if identifier != '':
            with open(self.data_path, 'r') as index:
                sequence = index.readline()
                while sequence:
                    if identifier in sequence:
                        param_seq = sequence.split(' ')
                        offset_head = int(param_seq[1])
                        offset = int(param_seq[2])
                        length = int(param_seq[3])
                        next_seq = index.readline()
                        while next_seq and identifier in next_seq:
                            length += int(next_seq.split(' ')[3])
                            next_seq = index.readline()
                        break
                    sequence = index.readline()
        return {'length': length, 'offset_head': offset_head, 'offset': offset}

    def get_sequences_of_range(self, min_range, max_range):
        sequences = []
        with open(self.data_path, 'r') as index:
            sequence = index.readline()
            while sequence and int(sequence.split(' ')[2]) < min_range:
                sequence = index.readline()

            while sequence and int(sequence.split(' ')[2]) < max_range:
                sequences.append(sequence)
                sequence = index.readline()
        return sequences



