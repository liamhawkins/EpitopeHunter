from operator import attrgetter

from tqdm import tqdm
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline


class Query:
    def __init__(self, filename):
        self.filename = filename
        seq_record = SeqIO.read(filename, 'fasta')
        self.id = seq_record.id
        self.seq = seq_record.seq

    def __len__(self):
        return len(self.seq)


class Subject:
    def __init__(self, filename, identity, buffer=0.01):
        self.filename = filename
        seq_record = SeqIO.read(filename, 'fasta')
        self.id = seq_record.id
        self.seq = seq_record.seq
        self.identity = identity
        self.identity_range = (identity-buffer, identity+buffer)

    def __str__(self):
        return self.id

    def in_range(self, best_hsp):
        if not best_hsp:
            return False
        identity = best_hsp.identities/len(best_hsp.query)
        return self.identity_range[0] <= identity <= self.identity_range[1]


class EpitopeFinder:
    def __init__(self, query, subjects, epitope_range=(10, 20), step_size=5):
        self.query = query
        self.subjects = subjects
        self.epitope_range = epitope_range
        self.step_size = step_size
        self.potential_epitopes = []
        self._run()

    def _run(self):
        for epitope_len in tqdm(range(self.epitope_range[0], self.epitope_range[1]+1)):
            self._scan_query(epitope_len)
        for e in self.potential_epitopes:
            print(e)

    def _scan_query(self, epitope_len):
        for start, end in self._gen_range(1, epitope_len, len(self.query), self.step_size):
            subjects_in_range = []
            for subject in self.subjects:
                best_hsp = self._blastp(query, (start, end), subject)
                subjects_in_range.append(subject.in_range(best_hsp))

            if all(subjects_in_range):
                self.potential_epitopes.append(f'{epitope_len} - {self.query.id}:{start}-{end}')

    @staticmethod
    def _blastp(query, query_loc, subject):
        NcbiblastpCommandline(query=query.filename, query_loc=f'{query_loc[0]}-{query_loc[1]}', subject=subject.filename, task='blastp-short', outfmt=5, out='seqs/temp.xml')()
        blast_record = NCBIXML.read(open('seqs/temp.xml'))
        try:
            return max(blast_record.alignments[0].hsps, key=attrgetter('score'))
        except IndexError:
            return None

    @staticmethod
    def _gen_range(start, length, stop, step):
        current_start = start
        current_end = start + length - 1
        while current_end <= stop:
            next_current_start = current_start + step
            next_current_end = current_end + step
            if next_current_end < stop:
                yield (current_start, current_end)
            else:
                yield (stop-length, stop)
            current_start = next_current_start
            current_end = next_current_end


if __name__ == '__main__':
    query = Query('seqs/epitope_range.fasta')
    subjects = [Subject('seqs/human_pfkfb1_1.fasta', 0.66),
                Subject('seqs/human_pfkfb4_1.fasta', 0.66),
                Subject('seqs/human_pfkfb4_2.fasta', 0.66),
                Subject('seqs/human_pfkfb3_1.fasta', 1.00)]
    EpitopeFinder(query, subjects)
