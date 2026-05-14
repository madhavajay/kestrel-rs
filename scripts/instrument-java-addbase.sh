#!/usr/bin/env bash
# Build an instrumented version of the Java kestrel.jar that logs per-addBase
# matrix bottom-row scores. Generates a `kestrel-instr.jar` next to the
# original `kestrel.jar` plus a per-J-R trace file.
#
# Usage (from kestrel-rs root):
#   scripts/instrument-java-addbase.sh OUTPUT_DIR
#
# Requires javac in PATH. The instrumentation modifies
# `kestrel/src/edu/gatech/kestrel/align/KmerAligner.java` in place — re-run
# `git checkout` to revert.
#
# The instrumented addBase emits a [JDBG-ADDBASE] trace line per call with
# fields: consensus_size, max_align_score, align_bot, gap_con_bot,
# max_pot_score, continue, base. Filter for J-R:4-119 via:
#
#   awk '/Building haplotypes.*J-R, start=4, end=119/{flag=1}
#        flag {print}
#        /Built.*\(fwd\): ActiveRegion\[name=J-R, start=4, end=119/{exit}' \
#       trace.log | grep JDBG-ADDBASE > jr-iter-trace.txt
#
# Then bisect against Rust's KDBG-CHOOSE trace from the J-R diagnostic test
# (run with `KESTREL_RUN_JR_DIAGNOSTIC=1 KESTREL_TRACE_REGION=J-R:4-119`).

set -euo pipefail

OUT_DIR=${1:-/tmp/kestrel-instrumented}
mkdir -p "$OUT_DIR"

KESTREL_ROOT=$(cd "$(dirname "$0")/.." && pwd)
SOURCE_FILE="$KESTREL_ROOT/kestrel/src/edu/gatech/kestrel/align/KmerAligner.java"
JAR_LIB="$KESTREL_ROOT/kestrel/lib"
CLASS_OUT="$OUT_DIR/classes"
mkdir -p "$CLASS_OUT"

if ! grep -q "JDBG-ADDBASE" "$SOURCE_FILE"; then
  echo "Patching KmerAligner.java with [JDBG-ADDBASE] instrumentation..."
  # Insert the trace before the return at the end of addBase.
  python3 - <<EOF
import re
path = "$SOURCE_FILE"
src = open(path).read()
needle = """tnSwap = matrixColGapCon;
\t\tmatrixColGapCon = matrixColGapConNext;
\t\tmatrixColGapConNext = tnSwap;
\t\t
\t\treturn maxPotScore >= maxAlignmentScore && maxPotScore > 0.0F;"""
insert = """tnSwap = matrixColGapCon;
\t\tmatrixColGapCon = matrixColGapConNext;
\t\tmatrixColGapConNext = tnSwap;

\t\tif (logger.isTraceEnabled()) {
\t\t\tfloat alignBot = (matrixColAlign[refLength - 1] != TraceNode.ZERO_NODE) ? matrixColAlign[refLength - 1].score : 0.0F;
\t\t\tfloat gapConBot = (matrixColGapCon[refLength - 1] != TraceNode.ZERO_NODE) ? matrixColGapCon[refLength - 1].score : 0.0F;
\t\t\tboolean cont = maxPotScore >= maxAlignmentScore && maxPotScore > 0.0F;
\t\t\tlogger.trace("[JDBG-ADDBASE] consensus_size={} max_align_score={} align_bot={} gap_con_bot={} max_pot_score={} continue={} base={}", consensusSize, maxAlignmentScore, alignBot, gapConBot, maxPotScore, cont, base);
\t\t}

\t\treturn maxPotScore >= maxAlignmentScore && maxPotScore > 0.0F;"""
if needle not in src:
    raise SystemExit("could not locate addBase return site to patch")
open(path, "w").write(src.replace(needle, insert))
print("Patched.")
EOF
fi

echo "Compiling..."
javac -cp "$JAR_LIB/kestrel.jar:$JAR_LIB/kanalyze.jar:$JAR_LIB/slf4j-api-1.7.12.jar" \
      -d "$CLASS_OUT" \
      "$SOURCE_FILE"

INSTR_JAR="$OUT_DIR/kestrel-instr.jar"
cp "$JAR_LIB/kestrel.jar" "$INSTR_JAR"
cd "$CLASS_OUT"
jar uf "$INSTR_JAR" \
    edu/gatech/kestrel/align/KmerAligner.class \
    edu/gatech/kestrel/align/KmerAligner\$AlignStart.class

echo "Instrumented jar at: $INSTR_JAR"
echo
echo "Run example (negative VNtyper FASTQ at caps 10/15):"
echo
echo "  java -cp '$INSTR_JAR:$JAR_LIB/kanalyze.jar:$JAR_LIB/slf4j-api-1.7.12.jar:$JAR_LIB/logback-core-1.1.3.jar:$JAR_LIB/logback-classic-1.1.3.jar:$JAR_LIB/java-getopt-1.0.14.jar:$JAR_LIB/commons-lang3-3.4.jar' \\"
echo "    edu.gatech.kestrel.clui.Main \\"
echo "    -k 20 --maxalignstates 10 --maxhapstates 15 \\"
echo "    -r REFERENCE_FASTA -o output.vcf -sNAME \\"
echo "    --temploc TMP_DIR --logstderr --loglevel TRACE \\"
echo "    R1.fastq R2.fastq 2> trace.log"
