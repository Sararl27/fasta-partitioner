import os
import unittest
import lithops
from pyfaidx import Fasta
import fastaPartitionerIndex as fp
import testsPartitionerFasta


def push_object_funct(local_input_path, data_bucket, input_data_prefix):
    bucket_objects = storage.list_keys(bucket=data_bucket)
    #storage.delete_objects(bucket=data_bucket, key_list=bucket_objects)
    for subdir, dirs, files in os.walk(local_input_path):
        print(subdir)
        for file_name in files:
            key = os.path.join(input_data_prefix, file_name)  # Added
            if key not in bucket_objects:  # Changed: if file_name not in bucket_objects:
                with open(os.path.join(subdir, file_name), 'rb') as file:  # Changed
                    print(f'\tUploading {key}...')
                    data = file.read()
                    storage.put_object(bucket=data_bucket, key=key, body=data)
                    print('\tOk!')
            else:  # Added
                print(f'\tIt is already uploaded: {key}...')  # Added
    print("Finished!")

    return storage.list_keys(bucket=data_bucket)


def generate_fasta_index_pyfaidx():
    Fasta('input_data/genes.fasta')
    print("Fasta index 'pyfaidx' done")


def generate_fasta_index_own(local_input_path, data_bucket, prefix, storage, key, workers):
    push_object_funct(local_input_path, data_bucket, prefix)
    fp.FastaPartitioner(storage, data_bucket, key, workers)


def test_partitioner_fasta():
    unittest.TextTestRunner(verbosity=2).run(unittest.TestLoader().loadTestsFromModule(testsPartitionerFasta))
    print(f'{testsPartitionerFasta.results[0]}:')
    for i in testsPartitionerFasta.results[1::]:
        print(f'\t{i}')


if __name__ == "__main__":
    storage = lithops.Storage()

    # Params
    data_bucket = 'fasta-partitioner'  # Change me
    prefix = 'fasta/'  # Change me
    local_input_path = './input_data'  # Change me

    key = f'fasta/genes.fasta'  # Change me
    obj = f'{data_bucket}/{key}'  # f'cos://{data_bucket}/{key}'  # Change me

    workers = 1000  # Change me

    # Execution
    # generate_fasta_index_pyfaidx()
    #generate_fasta_index_own(local_input_path, data_bucket, prefix, storage, key, workers)

    # run all tests with verbosity
    test_partitioner_fasta()
