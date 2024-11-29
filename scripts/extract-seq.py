"""Extracts regions to polish with additional flanks into fasta file"""
import argparse
import csv
import re
from collections import namedtuple
import btllib

Coordinate = namedtuple("Coordinate", "start end")
MIN_GAP_LENGTH = 1

def parse_args():
    """Parses arguments passed in command line"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--fasta", help="target file in fasta format", type=str, required=True
    )
    parser.add_argument(
        "-b",
        "--bed",
        help="bed file specifying regions to polish",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="prefix of output file [<output>.fa]",
        type=str,
        required=False,
        default="extracted_gaps",
    )
    parser.add_argument(
        "-l",
        "--length",
        help="length of flanking regions [64]",
        type=int,
        required=False,
        default=64,
    )
    return parser.parse_args()


def make_coord_dict(bed):
    """Creates a dictionary of tuples representing regions to polish"""
    coord_dict_2 = {}

    with open(bed, encoding="utf-8") as bed_file:
        bed_reader = csv.reader(bed_file, delimiter="\t", quotechar='"')
        for row in bed_reader:
            contig_name = row[0]
            coord = Coordinate(row[1], row[2])

            if contig_name not in coord_dict_2:
                coord_dict_2[contig_name] = [coord]
            else:
                coord_dict_2[contig_name].append(coord)
    return coord_dict_2

def extract_subsequences_from_bed(sequence, name, flank_length, writer, coords):
    """extracts seqs with coordinates from bed file and returns sequences with flank_length bp long flanks"""
    if name in coords:
        count = 0
        coord_list = coords[name]
        idx = 1

        filtered_coords = [] # short uppercase seqs appended to adj coordinates
        filtered_coords.append(
            coord_list[0]
        )  # first item always appended, avoids 0-index issues

        while idx < len(coord_list):
            coord = Coordinate(*coord_list[idx])
            prev_coord = Coordinate(*filtered_coords[-1])

            if (
                int(coord.start) - int(prev_coord.end)
            ) < 2 * flank_length:  # length between adjacent coords too small
                filtered_coords[-1] = (prev_coord.start, coord.end)
            else:
                filtered_coords.append(coord_list[idx])
            idx += 1

        for coord in filtered_coords:
            start = max(0, int(coord[0]) - flank_length)
            end = min(len(sequence), int(coord[1]) + flank_length)

            count += 1

            if start > end:
                write_flanked_subsequence(
                    sequence[start:end], start, end, name + "." + str(count), writer
                )


def write_flanked_subsequence(subsequence, start_flank, end_flank, gap_name, writer):
    """Writes subsequence with flanking regions into fasta file"""
    writer.write(gap_name, str(start_flank) + "-" + str(end_flank), subsequence.upper())


def main():
    """Parses fasta file to extract sequences + flanks"""
    args = parse_args()
    writer_fasta = btllib.SeqWriter(args.output, btllib.SeqWriter.FASTA)

    # makes coordinate dictionary if bed file provided
    coord_dict = make_coord_dict(args.bed)

    # loop through sequences in fasta file
    with btllib.SeqReader(args.fasta, btllib.SeqReaderFlag.LONG_MODE) as reader:
        for record in reader:
            seq_name, seq = record.id, str(record.seq)
            if seq_name in coord_dict:
                extract_subsequences_from_bed(
                    seq, seq_name, int(args.length), writer_fasta, coord_dict
                )

    writer_fasta.close()

if __name__ == "__main__":
    main()
