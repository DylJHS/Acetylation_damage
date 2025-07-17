#!/usr/bin/env python3
import sys
import os

if len(sys.argv) != 3:
    print("Usage: {} blast_output.txt total_sample_reads".format(sys.argv[0]))
    sys.exit(1)

blast_file = sys.argv[1]

try:
    total_sample_reads = float(sys.argv[2])
except ValueError:
    print("Error: total_sample_reads must be a number")
    sys.exit(1)

unique_queries = set()
total_hits = 0
sum_identity = 0.0
sum_mismatches = 0.0
sum_gap_openings = 0.0
sum_bit_score = 0.0

with open(blast_file, "r") as f:
    for line in f:
        if not line.strip():
            continue
        fields = line.strip().split("\t")
        if len(fields) < 12:
            continue  # Skip lines that don't have all expected columns
        query_id = fields[0]
        try:
            identity = float(fields[2])
            mismatches = float(fields[4])
            gap_openings = float(fields[5])
            bit_score = float(fields[11])
        except ValueError:
            continue  # Skip lines with non-numeric data in numeric columns
        unique_queries.add(query_id)
        total_hits += 1
        sum_identity += identity
        sum_mismatches += mismatches
        sum_gap_openings += gap_openings
        sum_bit_score += bit_score

alignment_rate = (len(unique_queries) / total_sample_reads) * 100 if total_sample_reads else 0.0
avg_identity = (sum_identity / total_hits) if total_hits else 0.0
avg_mismatches = (sum_mismatches / total_hits) if total_hits else 0.0
avg_gap_openings = (sum_gap_openings / total_hits) if total_hits else 0.0
avg_bit_score = (sum_bit_score / total_hits) if total_hits else 0.0

# Output the custom result file
output = (
    f"Alignment Rate: {alignment_rate:.2f}% ({len(unique_queries)} out of {int(total_sample_reads)} reads)\n"
    f"Total Number of Hits: {total_hits}\n"
    f"Average % Identity: {avg_identity:.2f}\n"
    f"Average Mismatches: {avg_mismatches:.2f}\n"
    f"Average Gap Openings: {avg_gap_openings:.2f}\n"
    f"Average Bit Score: {avg_bit_score:.2f}\n"
)

print(output)