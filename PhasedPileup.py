from pysam import Samfile
from collections import defaultdict
from os import path
import svgwrite as svg
from argparse import ArgumentParser
import pickle as pkl

def get_phase(read, snps):
    phase = None
    for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
        if ref_pos + 1 in snps:
            if phase == None:
                try:
                    # 1 if alternate, -1 if reference
                    phase = -1 + 2*snps[ref_pos + 1].index(read.seq[read_pos])
                except ValueError:
                    return 0 # This SNP isn't in the dataset
            else:
                try:
                    new_phase = -1 + 2*snps[ref_pos + 1].index(read.seq[read_pos])
                except ValueError:
                    return 0
                if new_phase != phase:
                    return 0 # read seems misphased
    return phase



def get_snps(snpfile):
    snps = defaultdict(dict)
    if path.exists('true_hets.tsv'):
        true_hets = {tuple(line.strip().split()):True
                     for line in open('true_hets.tsv')
                    }
    else:
        true_hets = defaultdict(lambda x: True)
    if path.exists(snpfile+'.pkl'):
        return pkl.load(open(snpfile+'.pkl', 'rb'))

    if snpfile.endswith('.bed'):
        for line in open(snpfile):
            chrom, _, start, refalt = line.strip().split()
            if true_hets.get((chrom, start), True):
                snps[chrom][int(start)] = refalt.split('|')
    else:
        raise NotImplementedError("We can't yet handle other SNP filetypes")

    pkl.dump(snps, open(snpfile+'.pkl', 'wb'))
    return snps

def get_snps_low_mem(snpfile, target_chrom, low, high):
    snps = defaultdict(dict)
    for line in open(snpfile):
        chrom, _, start, refalt = line.strip().split()
        if chrom  != target_chrom:
            continue
        start = int(start)
        if low <= start <= high:
            snps[chrom][start] = refalt.split('|')
    return snps

read_height = 8
x_scale = 1
spacing = 5
pos_color = 'blue'
neg_color = 'red'
unk_color = 'gray'
max_rows = 100


def insert_reads(read1, read2, reads_by_phase):
    if read2 is not None and read2.pos < read1.pos:
        read1, read2 = read2, read1
    for i, row in enumerate(reads_by_phase):
        if read1.reference_start > row[-1].reference_end + spacing:
            if read2 is None:
                row.append(read1)
            elif read1.pos < read2.pos:
                row.append(read1)
                row.append(read2)
            else:
                row.append(read2)
                row.append(read1)

            break
        elif i > max_rows:
            row.append(read1)
            if read2 is not None:
                row.append(read2)
            break
    else:
        reads_by_phase.append([read1])
        if read2 is not None:
            reads_by_phase[-1].append(read2)

def draw_read(read, dwg, x_start_coord, y_coord, phase_color, snps=None,
              with_snps_only=False, last_read=None):
    blocks = read.blocks
    slice_start, slice_end = 0, len(blocks)
    snp_locs = []
    for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
        if snps is not None and (ref_pos+1) in snps:
            snp_locs.append((read_pos, ref_pos))
    if with_snps_only and len(snp_locs) == 0:
        return
    g = dwg.g()
    g.add(svg.base.Title(("<tspan>QName {}: </tspan>"
                         " <tspan> IsRead1: {}</tspan>"
                          "<tspan>Blocks :{} </tspan>"
                          "<tspan>MPos: {} </tspan>"
                          "<tspan>:LastBlocks {}</tspan>"
                         ).format(
                             read.qname,
                             read.is_read1,
                             read.blocks,
                             read.mpos if read.is_proper_pair else "",
                             ("({})".format(last_read.blocks)
                              if ((last_read is not None) and (last_read.qname == read.qname))
                              else "")
                         )
                        )
         )
    if read.is_reverse:
        slice_start += 1
        block_start, block_end = blocks[0]
        last_stop = block_end
        g.add(dwg.polygon(
            [(x_scale * (block_start - x_start_coord), y_coord + read_height/2),
             (x_scale * (min(block_start+read_height/2, block_end) - x_start_coord), y_coord+read_height),
             (x_scale * (block_end - x_start_coord), y_coord+read_height),
             (x_scale * (block_end - x_start_coord), y_coord),
             (x_scale * (min(block_start+read_height/2, block_end) - x_start_coord), y_coord),
             (x_scale * (block_start - x_start_coord), y_coord+read_height/2)
            ],
            class_='read{}'.format(phase_color),
            id=read.qname,
        ))
    else:
        slice_end -= 1
        block_start, block_end = blocks[-1]
        g.add(dwg.polygon(
            [
                (x_scale * (block_start - x_start_coord), y_coord + read_height),
                (x_scale * (max(block_start, block_end - read_height/2) - x_start_coord), y_coord+read_height),
                (x_scale * (block_end - x_start_coord), y_coord+read_height/2),
                (x_scale * (max(block_start, block_end - read_height/2) - x_start_coord), y_coord),
                (x_scale * (block_start - x_start_coord), y_coord),
                (x_scale * (block_start - x_start_coord), y_coord+read_height),
            ],
            class_='read{}'.format(phase_color),
        ))
        last_stop = block_start

    for i in range(slice_start, slice_end):
        block_start, block_stop = blocks[i]
        g.add(dwg.line((x_scale*(last_stop - x_start_coord),
                          y_coord + 0.5 * read_height),
                         (x_scale*(block_start-x_start_coord),
                          y_coord + 0.5 * read_height),
                       class_='bar{}'.format(phase_color),
                         id=read.qname,
                        ))
        g.add(dwg.rect((x_scale*(block_start - x_start_coord),
                          y_coord),
                         (x_scale*(block_stop - block_start),
                          read_height),
                       class_='read{}'.format(phase_color),
                         id=read.qname,
                        ))
        last_stop = block_stop
    if snps:
        for read_pos, ref_pos in snp_locs:
            alleles = snps[ref_pos + 1]
            if read.seq[read_pos] in alleles:
                snp_color = [neg_color, pos_color][alleles.index(read.seq[read_pos])]
                snp_phase = ['neg', 'pos'][alleles.index(read.seq[read_pos])]
                if snp_phase == phase_color:
                    snp_color = 'black'
                    snp_phase = 'agree'
            else:
                snp_phase = 'unk'
                snp_color = 'black'
            g.add(dwg.line(
                (x_scale * (ref_pos - x_start_coord), y_coord),
                (x_scale * (ref_pos - x_start_coord), y_coord + read_height),
                class_='snp snp{}'.format(snp_phase),
            ))

    if last_read is not None and read.qname == last_read.qname:
        g.add(dwg.line((x_scale * (last_read.reference_end - x_start_coord),
                          y_coord + read_height/2),
                         (x_scale * (read.reference_start - x_start_coord),
                          y_coord + read_height/2),
                       class_='barinsert',
                        ))
    if read.is_read1:
        if read.is_reverse:
            read_start = blocks[-1][1]
            read_start_offset = read_start - 5
        else:
            read_start = blocks[0][0]
            read_start_offset = read_start + 5
        g.add(dwg.polygon(
            [
                (x_scale * (read_start - x_start_coord), y_coord + read_height),
                (x_scale * (read_start_offset - x_start_coord), y_coord+read_height/2),
                (x_scale * (read_start - x_start_coord), y_coord),
                (x_scale * (read_start - x_start_coord), y_coord+read_height),
            ],
            class_='r1_flag',
        ))
    dwg.add(g)


def find_transcript(annotation, splitkey=' '):
    annot_dict = dict(entry.strip().split(splitkey)
            for entry in annotation.strip(';"\'').split(';')
            if len(entry.strip().split()))
    return annot_dict.get('transcript_id', '???')


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('snp_file')
    parser.add_argument('samfile', nargs='+')
    parser.add_argument('--gene-name', '-g', default=None)
    parser.add_argument('--gtf-file', '-G', default=None)
    parser.add_argument('--coords', '-c', default=None,
                        help="Format: chrom:10..100 or chrom:10-100")
    parser.add_argument('--chrom', '-C', default=None)
    parser.add_argument('--draw-exons', '-x', default=False, action='store_true')
    parser.add_argument('--low-mem', default=True, action='store_true')
    parser.add_argument('--no-low-mem', dest='low_mem', action='store_false')
    parser.add_argument('--skip-if-deep', '-d', default=False,
                        action='store_true')
    parser.add_argument('--draw-all-snps', '-S', default=False,
                        action='store_true')
    parser.add_argument('--outdir', '-o', default='')
    args = parser.parse_args()

    if args.coords is not None:
        args.chrom, coords = args.coords.split(':')
        sep = '-' if '-' in coords else '..'
        coords = coords.replace(',', '').split(sep)

        args.coords = int(coords[0]), int(coords[-1])
    elif args.gene_name is not None:
        min_coord = 1e99
        max_coord = -1
        if args.gtf_file is None:
            args.gtf_file = 'mel_good.gtf'

        exons = set()
        for row in open(args.gtf_file):
            row = row.split('\t')
            if (('"{}"'.format(args.gene_name) in row[-1])
                    or ('={};'.format(args.gene_name) in row[-1])):
                min_coord = min(min_coord, int(row[3]))
                max_coord = max(max_coord, int(row[4]))
                args.chrom = row[0]
                curr_transcript = find_transcript(row[-1], '=' if
                        args.gtf_file.endswith('.gff') else ' ')
                exons.add((int(row[3]), int(row[4]), curr_transcript ))

        if max_coord <= min_coord:
            raise ValueError('Could not find gene: {} '.format(args.gene_name))
            assert False
        args.draw_exons = exons
        args.coords = min_coord, max_coord
    return args


def draw_exons(dwg, exons, x_start, y_start):
    transcript_dict = {}
    num_transcripts = 0
    for left, right, id in exons:
        if id not in transcript_dict:
            transcript_dict[id] = num_transcripts
            num_transcripts += 1
        curr_transcript = transcript_dict[id]
        #print(id, left, right, curr_transcript)

        dwg.add(dwg.rect(
            (x_scale * (left - x_start), y_start + read_height * curr_transcript),
            (x_scale * (right - left), read_height),
            style='opacity:1; fill: orange;',
                ))
        g = dwg.g(class_='hover_group')
        g.add(dwg.line((x_scale * (left-x_start), 0),
                       (x_scale * (left-x_start),
                        y_start + read_height * curr_transcript)))
        g.add(dwg.line((x_scale * (right-x_start), 0),
                       (x_scale * (right - x_start),
                        y_start + read_height * curr_transcript)))
        g.add(dwg.rect(
            (x_scale * (left - x_start), y_start + read_height * curr_transcript),
            (x_scale * (right - left), read_height),
            style='opacity:1; fill: orange;',
            ))
        y_start += read_height
        dwg.add(g)


if __name__ == "__main__":
    args = parse_args()
    phase_pos = []
    phase_neg = []
    phase_unk  = []

    # Note that phases are 0, 1, and -1
    phase_all = [phase_unk, phase_pos, phase_neg]
    if args.low_mem:
        snps = get_snps_low_mem(args.snp_file, args.chrom, *args.coords)
    else:
        snps = get_snps(args.snp_file)

    start_coord = 1e99
    end_coord = 0

    unmatched_reads = [{}, {}]
    # Note that in order to keep track of the un-phased reads, we need to do a
    # 2-pass approach to know the height of all of the classes
    gene_chrom = args.chrom
    gene_coords = args.coords
    num_reads = 0
    r1_lr = [0, 0]
    for fname in args.samfile:
        for read in (Samfile(fname).fetch(gene_chrom,
                                                 gene_coords[0],
                                                 gene_coords[1])):
            #if (not ((gene_coords[0] <= read.reference_start <= gene_coords[1])
                     #and (gene_coords[0] <= read.reference_end <= gene_coords[1]))):
                #continue
            for _, pos in read.get_aligned_pairs(matches_only=True):
                if gene_coords[0] <= pos <= gene_coords[1]: break
            else:
                continue
            phase = get_phase(read, snps[gene_chrom])
            num_reads += 1

            if read.is_read1:
                r1_lr[read.is_reverse] += 1

            # read.is_read1 = 1 if read_1, 0 if read_2
            if read.qname in unmatched_reads[read.is_read1]:
                other_read, other_phase = unmatched_reads[read.is_read1].pop(read.qname)

                if phase is None and other_phase is None:
                    phase = 0
                elif other_phase is None:
                    pass
                elif phase is None:
                    phase = other_phase
                elif phase != other_phase:
                    phase = 0
                else:
                    pass

                reads_by_phase = phase_all[phase]
                insert_reads(read, other_read, reads_by_phase)
            else:
                unmatched_reads[read.is_read2][read.qname] = (
                    read, phase
                )


            if phase is None:
                continue

            start_coord = min(start_coord, read.reference_start)
            end_coord = max(end_coord, read.reference_end)

    for reads in unmatched_reads:
        for qname in reads:
            read, phase = reads[qname]
            if phase is None:
                continue
            reads_by_phase = phase_all[phase]
            insert_reads(read, None, reads_by_phase)


    print("Matched ", num_reads)


    max_depth_pos = len(phase_pos)
    max_depth_neg = len(phase_neg)
    max_depth_unk = len(phase_unk)

    out_fname = args.samfile[0].replace('.bam', '_phased.svg')
    if args.gene_name:
        out_fname = path.join(path.dirname(args.outdir or out_fname),
                args.gene_name + '_' + path.basename(out_fname))
    start_coord -= 10
    end_coord += 10
    if args.draw_exons:
        start_coord = min(gene_coords[0] - 10, start_coord)
        end_coord = max(gene_coords[1] + 10, end_coord)
    dwg = svg.Drawing(out_fname,
                      (x_scale*(end_coord - start_coord),
                       read_height * 1.2 * (max_depth_pos +
                                            max_depth_neg +
                                            len(args.draw_exons) +
                                            max_depth_unk +  15)),
                      profile='full',
                      )
    dwg.add(dwg.style(
        (
            '.hover_group{{opacity:0;}} \n'
            '.hover_group:hover \n{{\n\topacity:1;\n\tstroke-width:1!important;\n\tstroke:#000000;\n}}'
            '.barneg{{stroke-width:1;stroke:{neg};}}\n'
            '.barpos{{stroke-width:1;stroke:{pos};}}\n'
            '.barunk{{stroke-width:1;stroke:{unk};}}\n'
            '.barinsert{{stroke-width:1;stroke:black;}}\n'
            '.readneg{{fill:{neg};}}\n'
            '.readpos{{fill:{pos};}}\n'
            '.readunk{{fill:{unk};}}\n'
            '.snpagree{{stroke:black}}\n'
            '.snpneg{{stroke:{neg};}}\n'
            '.snppos{{stroke:{pos};}}\n'
            '.snpunk{{stroke:white}}\n'
            '.snp{{stroke-width:1}}\n'
            '.r1flag{{fill:black;}}\n'
        ).format(**{'pos': pos_color, 'neg': neg_color, 'unk': unk_color})
    ))

    y_start = 10 + max_depth_neg * 1.2*read_height

    last_read = None
    for row_num, row in enumerate(phase_neg):
        if args.skip_if_deep and row_num == (max_rows+1):
            #print("bailing")
            continue
        for read in row:
            if read is None:
                continue
            draw_read(read, dwg,
                      start_coord, y_start - row_num * 1.2 * read_height,
                      'neg',
                      snps=snps[gene_chrom],
                      last_read = last_read
                     )
            last_read = read

    y_start += 3 * read_height
    for row_num, row in enumerate(phase_unk):
        if args.skip_if_deep and row_num == (max_rows+1):
            continue
        for read in row:
            draw_read(read, dwg,
                      start_coord, y_start + row_num * 1.2 * read_height,
                      'unk',
                      snps=snps[gene_chrom],
                      last_read = last_read
                     )
            last_read = read
    y_start += max_depth_unk * 1.2 * read_height + 3 * read_height

    for row_num, row in enumerate(phase_pos):
        if args.skip_if_deep and row_num == (max_rows+1):
            continue
        for read in row:
            draw_read(read, dwg,
                      start_coord, y_start + row_num * 1.2 * read_height,
                      'pos',
                      snps=snps[gene_chrom],
                      last_read = last_read
                     )
            last_read = read
    y_start += max_depth_pos * 1.2 * read_height + 3 * read_height

    if args.draw_all_snps:
        g = dwg.g()
        for snp in snps[gene_chrom]:
            if snp < start_coord or snp > end_coord: continue
            g.add(dwg.line(
                (x_scale * (snp - start_coord), 10 ),
                (x_scale * (snp - start_coord), y_start),
                class_='snp snp{}'.format('agree'),
            ))
        dwg.add(g)



    if args.draw_exons:
        draw_exons(dwg, args.draw_exons, start_coord, y_start)
    dwg.save()


