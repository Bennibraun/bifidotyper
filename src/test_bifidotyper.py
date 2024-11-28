import unittest
import os
from unittest.mock import patch
from bifidotyper import parse_arguments, encode_fastq_files

class TestBifidotyper(unittest.TestCase):

    def test_parse_arguments_single_end(self):
        """Test single-end FASTQ file argument parsing."""
        test_args = ['bifidotyper.py', '--single-end', 'sample1.fastq', 'sample2.fastq']
        with patch('sys.argv', test_args):
            args = parse_arguments()
            self.assertEqual(args.single_end, ['sample1.fastq', 'sample2.fastq'])
            self.assertIsNone(args.paired_end)

    def test_parse_arguments_paired_end(self):
        """Test paired-end FASTQ file argument parsing."""
        test_args = [
            'bifidotyper.py', 
            '--paired-end', 'read1.fastq', 'read2.fastq',
            '--paired-end', 'read3.fastq', 'read4.fastq'
        ]
        with patch('sys.argv', test_args):
            args = parse_arguments()
            self.assertIsNone(args.single_end)
            self.assertEqual(
                args.paired_end, 
                [['read1.fastq', 'read2.fastq'], ['read3.fastq', 'read4.fastq']]
            )

    def test_parse_arguments_both_single_and_paired(self):
        """Test both single-end and paired-end arguments."""
        test_args = [
            'bifidotyper.py',
            '--single-end', 'sample1.fastq', 'sample2.fastq',
            '--paired-end', 'read1.fastq', 'read2.fastq'
        ]
        with patch('sys.argv', test_args):
            args = parse_arguments()
            self.assertEqual(args.single_end, ['sample1.fastq', 'sample2.fastq'])
            self.assertEqual(args.paired_end, [['read1.fastq', 'read2.fastq']])

    def test_encode_fastq_files_single_end(self):
        """Test encoding of single-end FASTQ files."""
        single_end_files = ['sample1.fastq', 'sample2.fastq']
        expected = {
            'sample1.fastq': os.path.abspath('sample1.fastq'),
            'sample2.fastq': os.path.abspath('sample2.fastq'),
        }
        result = encode_fastq_files(single_end_files, paired_end_files=None)
        self.assertIn('single_end', result)
        self.assertEqual(result['single_end'], expected)

    def test_encode_fastq_files_paired_end(self):
        """Test encoding of paired-end FASTQ files."""
        paired_end_files = [
            ('read1.fastq', 'read2.fastq'),
            ('read3.fastq', 'read4.fastq')
        ]
        expected = {
            'read1.fastq': (os.path.abspath('read1.fastq'), os.path.abspath('read2.fastq')),
            'read3.fastq': (os.path.abspath('read3.fastq'), os.path.abspath('read4.fastq')),
        }
        result = encode_fastq_files(single_end_files=None, paired_end_files=paired_end_files)
        self.assertIn('paired_end', result)
        self.assertEqual(result['paired_end'], expected)

    def test_encode_fastq_files_both_single_and_paired(self):
        """Test encoding of both single-end and paired-end FASTQ files."""
        single_end_files = ['sample1.fastq']
        paired_end_files = [('read1.fastq', 'read2.fastq')]
        expected = {
            'single_end': {
                'sample1.fastq': os.path.abspath('sample1.fastq'),
            },
            'paired_end': {
                'read1.fastq': (os.path.abspath('read1.fastq'), os.path.abspath('read2.fastq')),
            }
        }
        result = encode_fastq_files(single_end_files, paired_end_files)
        self.assertEqual(result, expected)

if __name__ == '__main__':
    unittest.main()
