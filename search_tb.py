from six_frames import generate_six_frames
import re,os
from Bio import SeqIO
hmmer_path = '/usr/bin/hmmsearch'

def generate_fasta(six_frames):
    all_str = ''
    for seqid,offset,seq in six_frames:
        new_seqid = str(offset)
        all_str += '>%s\n%s\n' % (new_seqid,seq)
    return all_str

def convert_hmmer_format(tmp):
    start_line = 'Domain annotation for each sequence:'
    end_line = "Internal pipeline statistics summary:"
    tmp = tmp[tmp.index(start_line):tmp.index(end_line)]

    info_dict = {}
    all_lines = tmp.split('\n>>')
    for anno in all_lines:
        anno = anno.split('\n')
        if not anno:
            continue
        anno = [_ for _ in anno if _]
        if (not anno[0].startswith(start_line)) and (not anno[0].startswith(end_line)):
            name = anno[0].strip(' ').split(' ')[0]
            header = anno[1].strip(' ').split(' ')
            header = [_ for _ in header if _]
            values = []
            for _ in anno[3:]:
                value = _.strip(' ').split(' ')
                value = [_ for _ in value if _]
                value = (float(value[5]),float(value[-3]), float(value[-4]))
                values.append(value)

            info_dict[name] = list(sorted(values))[0][1:3]
    return info_dict

def exec_hmmer(six_frames,hmmer_profile):
    tmp = generate_fasta(six_frames)
    with open('/home/liaoth/tmp_180712','w') as f1:
        f1.write(tmp)
    cmdline = "%s --noali '%s' '%s'" % (hmmer_path,hmmer_profile,'/home/liaoth/tmp_180712')
    result = os.popen(cmdline)
    result = result.read()
    # print(cmdline)
    os.system('rm /home/liaoth/tmp_180712')
    return result

def parse2real_pos(result,ori_fasta):
    if len(result) != 1:
        raise IOError
    fasta = list(SeqIO.parse(ori_fasta,format='fasta'))[0]

    key = list(result.keys())[0]
    if 'rev' in key:
        seq = str(fasta.reverse_complement().seq).upper()
    else:
        seq = str(fasta.seq).upper()
    offset = int(key[-1])
    nucl_pos = [(_ + 1)*3 + offset -1  for _ in result[key]]

    return nucl_pos,seq[int(min(nucl_pos)):int(max(nucl_pos)+1)]


if __name__ == '__main__':
    import glob
    six_frames = generate_six_frames('/home/liaoth/data2/project/NR-Pfam/searching_all_species/updated_180525/cluster_all_nasN/raw_data/nucl/CP000926.1.fasta')
    hmmer_profile = "/home/liaoth/data2/project/NR-Pfam/searching_all_species/updated_180525/cluster_all_nasN/mafft_result(protein)/domains_hmm/Molydop_binding.hmm"
    tmp = exec_hmmer(six_frames,hmmer_profile)
    result = convert_hmmer_format(tmp)

    ############################################################
    # extract linker
    f1 = open('/home/liaoth/data2/project/NR-Pfam/searching_all_species/updated_180525/cluster_all_nasN/raw_data/from_pro2nucl/extract_linker.fasta','w')
    for i in glob.glob('/home/liaoth/data2/project/NR-Pfam/searching_all_species/updated_180525/cluster_all_nasN/raw_data/from_pro2nucl/pro2nucl/*.fasta'):
        if 'nasN' not in os.path.basename(i):
            result = convert_hmmer_format(exec_hmmer(generate_six_frames(i), hmmer_profile))
            result2 = convert_hmmer_format(exec_hmmer(generate_six_frames(i), '/home/liaoth/data2/project/NR-Pfam/searching_all_species/updated_180525/cluster_all_nasN/mafft_result(protein)/domains_hmm/Flavodoxin_1.hmm'))
            if list(result.keys())[0] == list(result2.keys())[0]:
                key = list(result.keys())[0]
                pos = max(result[key]),min(result2[key])
                seq_extract = parse2real_pos({key:pos},i)[1]
                # tmp = list(SeqIO.parse(i,format='fasta'))[0].description
                tmp = os.path.basename(i).split('.fasta')[0]
                f1.write('>%s\n%s\n' % (total_nasN_dict[tmp],seq_extract))
                f1.flush()
            else:
                print('error')
    ############################################################
    # extract nasA,nasN
    f1 = open('/home/liaoth/data2/project/NR-Pfam/searching_all_species/updated_180525/cluster_all_nasN/raw_data/from_pro2nucl/extract_nasA.fasta','w')
    f2 = open('/home/liaoth/data2/project/NR-Pfam/searching_all_species/updated_180525/cluster_all_nasN/raw_data/from_pro2nucl/extract_nasN.fasta','w')
    for i in glob.glob('/home/liaoth/data2/project/NR-Pfam/searching_all_species/updated_180525/cluster_all_nasN/raw_data/from_pro2nucl/pro2nucl/*.fasta'):
        if 'nasN' not in os.path.basename(i):
            result = convert_hmmer_format(exec_hmmer(generate_six_frames(i), hmmer_profile))
            result2 = convert_hmmer_format(exec_hmmer(generate_six_frames(i), '/home/liaoth/data2/project/NR-Pfam/searching_all_species/updated_180525/cluster_all_nasN/mafft_result(protein)/domains_hmm/Flavodoxin_1.hmm'))
            if list(result.keys())[0] == list(result2.keys())[0]:
                key = list(result.keys())[0]
                pos = max(result[key]),0
                seq_extract = parse2real_pos({key:pos},i)[1]
                # tmp = list(SeqIO.parse(i,format='fasta'))[0].description
                tmp = os.path.basename(i).split('.fasta')[0]
                f1.write('>%s\n%s\n' % (total_nasN_dict[tmp],seq_extract))
                f1.flush()
                ############################################################
                pos = min(result2[key]),len(list(SeqIO.parse(i,format='fasta'))[0])-1
                seq_extract = parse2real_pos({key:pos},i)[1]
                # tmp = list(SeqIO.parse(i,format='fasta'))[0].description
                tmp = os.path.basename(i).split('.fasta')[0]
                f2.write('>%s\n%s\n' % (total_nasN_dict[tmp],seq_extract))
                f2.flush()
            else:
                print('error')
