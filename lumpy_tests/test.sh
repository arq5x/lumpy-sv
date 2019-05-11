test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

set -o nounset


SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

export LUMPY_HOME=$DIR/../lumpy-sv/
export LUMPY=$DIR/../bin/lumpy

cd $DIR/orig

wget -nc https://s3.amazonaws.com/lumpy/pe.pos_sorted.bam
wget -nc https://s3.amazonaws.com/lumpy/sr.pos_sorted.bam
wget -nc https://s3.amazonaws.com/lumpy/sim.bedpe

run original_pe_test ./bp_pe.sh 1 20 1e-3 150 pe.pos_sorted.bam 
assert_in_stderr chr10

if [ `command -v bedtools` ];then
    assert_equal 40 $(cat $STDOUT_FILE | wc -l)

    assert_equal 40 $( bedtools pairtopair \
                        -a sim.bedpe \
                        -b $STDOUT_FILE \
                        -type both -is -slop 3 \
                        | cut -f7 | sort -u | wc -l)
    assert_equal 0 $( bedtools pairtopair \
                        -b sim.bedpe \
                        -a $STDOUT_FILE \
                        -type notboth -is -slop 3 \
                        | cut -f7 | sort -u | wc -l)
fi

run original_sr_test ./bp_sr.sh 1 20 1e-3 150 sr.pos_sorted.bam
assert_in_stderr chr10

if [ `command -v bedtools` ];then
    assert_equal 44 $(cat $STDOUT_FILE | wc -l)

    assert_equal 43 $( bedtools pairtopair \
                        -a sim.bedpe \
                        -b $STDOUT_FILE \
                        -type both -is -slop 3 \
                        | cut -f7 | sort -u | wc -l)
    assert_equal 1 $( bedtools pairtopair \
                        -b sim.bedpe \
                        -a $STDOUT_FILE \
                        -type notboth -is -slop 3 \
                        | cut -f7 | sort -u | wc -l)
fi

run original_pesr_test ./bp_pesr.sh 1 20 1e-3 150 pe.pos_sorted.bam sr.pos_sorted.bam 
assert_in_stderr chr10

if [ `command -v bedtools` ];then
    assert_equal 276 $(cat $STDOUT_FILE | wc -l)

    assert_equal 275 $( bedtools pairtopair \
                        -a sim.bedpe \
                        -b $STDOUT_FILE \
                        -type both -is -slop 3 \
                        | cut -f7 | sort -u | wc -l)
    assert_equal 1 $( bedtools pairtopair \
                        -b sim.bedpe \
                        -a $STDOUT_FILE \
                        -type notboth -is -slop 3 \
                        | cut -f7 | sort -u | wc -l)
fi

cd ..

cd rice 

wget -nc https://s3.amazonaws.com/lumpy/AL87.discordant.sort.bam
wget -nc https://s3.amazonaws.com/lumpy/AL87.histo
wget -nc https://s3.amazonaws.com/lumpy/AL87.sr.sort.bam
wget -nc https://s3.amazonaws.com/lumpy/rice.pesr.vcf

run rice bash run.sh
assert_in_stderr Chr1
assert_in_stderr ChrUn
assert_equal 1550 $(cat $STDOUT_FILE | wc -l)
assert_equal \
    0 \
    $(bedtools intersect -v -a rice.pesr.vcf -b $STDOUT_FILE | wc -l)

cd ..

cd honeybee 

wget -nc https://s3.amazonaws.com/lumpy/test.filtered_split.bam
wget -nc https://s3.amazonaws.com/lumpy/test.sample1.lib1.x4.histo
wget -nc https://s3.amazonaws.com/lumpy/test_both.bam.vcf
wget -nc https://s3.amazonaws.com/lumpy/test_onlydisc.bam.vcf
wget -nc https://s3.amazonaws.com/lumpy/test_onlysplit.bam.vcf
wget -nc https://s3.amazonaws.com/lumpy/test_sub.filtered_disc.bam

run honeybee_pe bash run_pe.sh
assert_equal \
    0 \
    $(bedtools intersect -v -a test_onlydisc.bam.vcf -b $STDOUT_FILE | wc -l)

run honeybee_sr bash run_sr.sh
assert_equal \
    0 \
    $(bedtools intersect -v -a test_onlysplit.bam.vcf -b $STDOUT_FILE | wc -l)

run honeybee_pesr bash run_pesr.sh
assert_equal \
    0 \
    $(bedtools intersect -v -a test_both.bam.vcf -b $STDOUT_FILE | wc -l)

cd ..


cd NA12878

wget -nc https://s3.amazonaws.com/lumpy/NA12878.NA12891.NA12892.vcf.o
wget -nc https://s3.amazonaws.com/lumpy/NA12878.bam
wget -nc https://s3.amazonaws.com/lumpy/NA12878.bam.bai
wget -nc https://s3.amazonaws.com/lumpy/NA12891.bam
wget -nc https://s3.amazonaws.com/lumpy/NA12891.bam.bai
wget -nc https://s3.amazonaws.com/lumpy/NA12892.bam
wget -nc https://s3.amazonaws.com/lumpy/NA12892.bam.bai

rm -f ./*.{disc,split}.bam
run lumpy_smooth ../../scripts/lumpy_smooth NA12878.bam NA12891.bam NA12892.bam
assert_exit_code 0 \
assert_equal \
    5 \
    $( grep -v "^#" NA12878.NA12891.NA12892.vcf | cut -f 10- | wc -l )


assert_equal \
    0 \
    $(diff <(grep -v "^#" NA12878.NA12891.NA12892.vcf.o | cut -f 10- ) <(grep -v "^#" NA12878.NA12891.NA12892.vcf | cut -f 10- ) | wc -l)

run lumpy_smooth_reuse ../../scripts/lumpy_smooth -n trio NA12878.bam NA12891.bam NA12892.bam
assert_exit_code 0
assert_equal $(grep -c "using existing splitters:" $STDOUT_FILE) 3

assert_equal 5 $(zgrep -cv ^# trio.vcf)
assert_equal 5 $(zgrep -cv ^# trio.svtyped.vcf.gz)

cd ..


cd lumpyexpress

rm -f ./*.vcf
run lumpyexpress ../../scripts/lumpyexpress -B NA12878.bam -S NA12878.split.bam -D NA12878.disc.bam -o NA12878.vcf
assert_exit_code 0
assert_equal \
    0 \
    $(diff <(grep -v '^#' NA12878.vcf.o ) <(grep -v '^#' NA12878.vcf) | wc -l)

run lumpyexpress ../../scripts/lumpyexpress -B NA12878.w_space.bam -S NA12878.w_space.split.bam -D NA12878.w_space.disc.bam -o NA12878.w_space.vcf
assert_exit_code 0
assert_equal \
    0 \
    $(diff <(grep -v '^#' NA12878.vcf.o ) <(grep -v '^#' NA12878.w_space.vcf) | wc -l)

cd ..
