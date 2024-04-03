process racon_polish {
    publishDir 'tmp'
    label 'racon'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input: tuple val(id), path(draft_fastas), path(fastqs)
    output: tuple val(id), path("racond",type:'dir')
    shell: '''
mkdir -p pafs
find -L !{draft_fastas} -name "*.fasta" \
    | parallel -j !{task.cpus} \
        'minimap2 -t 1 {} !{fastqs}/{/.}.fastq > pafs/{/.}.paf'
mkdir -p racond
# Note that these racon parameters are as expected for medaka! Don't change them
find -L pafs -name "*.paf" -size +1b \
    | parallel -j !{task.cpus} ' \
        racon -m 8 -x -6 -g -8 -w 500 \
            -t 1 \
            !{fastqs}/{/.}.fastq \
            {} \
            !{draft_fastas}/{/.}.fasta \
            -u \
            > racond/{/.}.fasta \
        || echo "Error with run {/.}" \
        '
'''
}

process medaka_polish {
    cache 'lenient'
    publishDir 'tmp'
    label 'medaka'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input: tuple val(id), path(draft_fastas), path(fastqs), path(medaka_model)
    output: tuple val(id), path("medakd",type:'dir')
    shell: '''
export TF_FORCE_GPU_ALLOW_GROWTH=true
mkdir -p aligns
cp -Lr !{draft_fastas} copy_!{draft_fastas}
find -L copy_!{draft_fastas} -name "*.fasta" \
    | parallel \
        -j !{Math.round(task.cpus)} \
            'mini_align -m -r {} -i !{fastqs}/{/.}.fastq -p aligns/{/.} \
                || echo "Error with run {}" '
# -t threads
# get model
#wget "https://media.githubusercontent.com/media/nanoporetech/medaka/master/medaka/data/medaka_config_model.tar.gz"
# consensus
mkdir -p results
export TF_FORCE_GPU_ALLOW_GROWTH=true 
export TIMEOUT=300
find -L aligns -name "*.bam" \
    | parallel --progress --delay 0.2 --timeout $TIMEOUT --joblog consensus.log \
        -j !{Math.round(task.cpus/2)} \
            'exit 1 & medaka consensus {} results/{/.}.hdf \
            --chunk_len 100 --chunk_ovlp 10 \
            --model !{medaka_model} \
            --batch 1 --threads 2' \
    || echo "ok"
parallel --progress --joblog consensus.log --retry-failed --timeout $TIMEOUT \
    -j !{Math.round(task.cpus/4)}
parallel --progress --joblog consensus.log --retry-failed --timeout $TIMEOUT \
    -j !{Math.round(task.cpus/4)}
parallel --progress --joblog consensus.log --retry-failed --timeout $TIMEOUT \
    -j !{Math.round(task.cpus/8)}
parallel --progress --joblog consensus.log --retry-failed --timeout $TIMEOUT \
    -j !{Math.round(task.cpus/16)}
parallel --progress --joblog consensus.log --retry-failed --timeout $TIMEOUT \
    -j 1
# stitch
mkdir -p medakd
find -L results -name "*.hdf" \
    | parallel --progress --timeout $TIMEOUT --joblog stitch.log \
        -j !{Math.round(task.cpus)} \
            'exit 1 & medaka stitch {} copy_!{draft_fastas}/{/.}.fasta medakd/{/.}.fasta' \
    || echo "ok"
parallel --progress --joblog stitch.log --retry-failed --timeout $TIMEOUT \
    -j 1
sync # trying this out to see if it stops reruns...
sleep 10 # this is to try and prevent from rerunning due to weird caching???
# TODO disable see if it's still a problem
'''
}

