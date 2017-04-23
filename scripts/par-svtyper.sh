set -euo pipefail

export vcf=$1
export bams=$2

export LINES=250

svtpar() {
    set -euo pipefail
    subvcf=$1
    lib=$TMPDIR/tmp.$mypid.$(echo $bams | md5sum | awk '{ print $1 }').lib
    svtyper -B $bams -i $subvcf -l $lib
    # vcf is a temporary file so clean it up.
    # but for the first time, it's and fd from the header so the rm will fail.
    set +e
    tail -1 $subvcf 2>/dev/null | awk '{ print "svtyped through:" $1":"$2 }' > /dev/stderr
    rm -f $subvcf >/dev/null 2>&1
}
export -f svtpar

genvcfs() {
    set -euo pipefail
    lines=$1
    vcf=$2
    set +e
    header="$(less $vcf | head -2000 | grep ^#)"
    set -e
    # run this once to force the calculation of the bam stats.
    # this avoids a race condition.

    zless "$vcf" | awk -vtmp=$TMPDIR -vtmpvcf=xx -vpid=$$ -v header="$header"  -vLINES=$lines 'BEGIN{nth=0;} !(NR % LINES) {
        # we print the tmpvcf and move onto the next.
        if(NR > 1) {
            print tmpvcf;
        }
        tmpvcf=tmp"/tmp"pid"."nth".vcf"
        print header > tmpvcf;
        nth+=1;
    }
    {
        print $0 > tmpvcf;
    }
    END {
        print tmpvcf;
    }' 
}

vcfsort() {
     awk 'BEGIN{x=0;} $0 ~/^#/{ if(x==0) {print;} if($1 == "#CHROM"){x=1}; next;} {x=1; print $0 | "LC_ALL=C sort --buffer-size 2G -k1,1 -k2,2n" }'
}


export -f genvcfs

if [[ -z ${TMPDIR+x} ]]; then
       export TMPDIR=/tmp/
fi
export mypid=$$

cleanup() {
   rm -f $TMPDIR/tmp.$mypid.*.vcf
   rm -f $TMPDIR/tmp.$mypid.*.lib
}

trap cleanup EXIT

# add gargs to path
export PATH=$PATH:.

# get gargs if it doesn't exist
type gargs >/dev/null 2>&1 || {
    if [ "$(uname)" == "Darwin" ]; then
        wget -O gargs https://github.com/brentp/gargs/releases/download/v0.3.5/gargs_darwin
    else
        wget -O gargs https://github.com/brentp/gargs/releases/download/v0.3.5/gargs_linux
    fi
    chmod +x ./gargs
}

psvtyper() {
    set -euo pipefail
    export header=$(head -2000 $vcf | grep ^#)
    svtpar <(echo -e "$header")
    genvcfs $LINES $vcf | gargs -p $THREADS "svtpar {}" | vcfsort
}

mkdir -p $TMPDIR

export -f psvtyper
psvtyper
