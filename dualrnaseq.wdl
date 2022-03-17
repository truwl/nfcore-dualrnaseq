version 1.0

workflow dualrnaseq {
	input{
		String? name
		File samplesheet = "data/*{1,2}.fastq.gz"
		Boolean? single_end
		String outdir = "./results"
		Int max_cpus = 16
		String max_memory = "128.GB"
		String max_time = "240.h"
		String? fasta_host
		String? fasta_pathogen
		String? gff_host_genome
		String? gff_host_tRNA
		String? gff_pathogen
		String? transcriptome_host
		String? transcriptome_pathogen
		Boolean? read_transcriptome_fasta_host_from_file
		Boolean? read_transcriptome_fasta_pathogen_from_file
		String genome_host = "GRCh38"
		String genome_pathogen = "SL1344"
		Boolean? genomes_ignore
		Boolean? skip_fastqc
		String? fastqc_params
		Boolean? run_cutadapt
		String a = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
		String A = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
		Int quality_cutoff = 10
		String? cutadapt_params
		Boolean? run_bbduk
		Int minlen = 18
		String qtrim = "r"
		Int trimq = 10
		String ktrim = "r"
		Int k = 17
		Int mink = 11
		Int hdist = 1
		String adapters = "data/adapters.fa"
		String? bbduk_params
		String? libtype
		Int? incompatPrior
		Boolean? generate_salmon_uniq_ambig
		String gene_feature_gff_to_create_transcriptome_host = "['exon', 'tRNA']"
		String gene_feature_gff_to_create_transcriptome_pathogen = "['gene', 'sRNA', 'tRNA', 'rRNA']"
		String gene_attribute_gff_to_create_transcriptome_host = "transcript_id"
		String gene_attribute_gff_to_create_transcriptome_pathogen = "locus_tag"
		Boolean? run_salmon_selective_alignment
		Int kmer_length = 21
		Boolean? writeUnmappedNames
		Boolean? softclipOverhangs
		Boolean? dumpEq
		Boolean? writeMappings
		Boolean? keepDuplicates
		String? salmon_sa_params_index
		String? salmon_sa_params_mapping
		Boolean? run_salmon_alignment_based_mode
		String? salmon_alignment_based_params
		Boolean? run_star
		String outSAMunmapped = "Within"
		String outSAMattributes = "Standard"
		Int outFilterMultimapNmax = 999
		String outFilterType = "BySJout"
		Int alignSJoverhangMin = 8
		Int alignSJDBoverhangMin = 1
		Int outFilterMismatchNmax = 999
		Int outFilterMismatchNoverReadLmax = 1
		Int alignIntronMin = 20
		Int alignIntronMax = 1000000
		Int alignMatesGapMax = 1000000
		Int limitBAMsortRAM = 0
		Int winAnchorMultimapNmax = 999
		Int sjdbOverhang = 100
		String quantTranscriptomeBan = "Singleend"
		String? star_salmon_index_params
		String? star_salmon_alignment_params
		String outWigType = "None"
		String outWigStrand = "Stranded"
		String? star_index_params
		String? star_alignment_params
		Boolean? run_htseq_uniquely_mapped
		String stranded = "yes"
		Int max_reads_in_buffer = 30000000
		Int minaqual = 10
		String? htseq_params
		String gene_feature_gff_to_quantify_host = "['exon', 'tRNA']"
		String host_gff_attribute = "gene_id"
		String gene_feature_gff_to_quantify_pathogen = "['gene', 'sRNA', 'tRNA', 'rRNA']"
		String pathogen_gff_attribute = "locus_tag"
		Boolean? mapping_statistics
		String rna_classes_to_replace_host = "{base_dir}/data/RNA_classes_to_replace.csv"
		String publish_dir_mode = "copy"
		String? email_on_fail
		Boolean? plaintext_email
		String max_multiqc_email_size = "25.MB"
		Boolean? monochrome_logs
		String? multiqc_config
		String tracedir = "./results/pipeline_info"
		Boolean? help
		String? email
		String custom_config_version = "master"
		String custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/master"
		String? hostnames
		String? config_profile_description
		String? config_profile_contact
		String? config_profile_url

	}

	call make_uuid as mkuuid {}
	call touch_uuid as thuuid {
		input:
			outputbucket = mkuuid.uuid
	}
	call run_nfcoretask as nfcoretask {
		input:
			name = name,
			samplesheet = samplesheet,
			single_end = single_end,
			outdir = outdir,
			max_cpus = max_cpus,
			max_memory = max_memory,
			max_time = max_time,
			fasta_host = fasta_host,
			fasta_pathogen = fasta_pathogen,
			gff_host_genome = gff_host_genome,
			gff_host_tRNA = gff_host_tRNA,
			gff_pathogen = gff_pathogen,
			transcriptome_host = transcriptome_host,
			transcriptome_pathogen = transcriptome_pathogen,
			read_transcriptome_fasta_host_from_file = read_transcriptome_fasta_host_from_file,
			read_transcriptome_fasta_pathogen_from_file = read_transcriptome_fasta_pathogen_from_file,
			genome_host = genome_host,
			genome_pathogen = genome_pathogen,
			genomes_ignore = genomes_ignore,
			skip_fastqc = skip_fastqc,
			fastqc_params = fastqc_params,
			run_cutadapt = run_cutadapt,
			a = a,
			A = A,
			quality_cutoff = quality_cutoff,
			cutadapt_params = cutadapt_params,
			run_bbduk = run_bbduk,
			minlen = minlen,
			qtrim = qtrim,
			trimq = trimq,
			ktrim = ktrim,
			k = k,
			mink = mink,
			hdist = hdist,
			adapters = adapters,
			bbduk_params = bbduk_params,
			libtype = libtype,
			incompatPrior = incompatPrior,
			generate_salmon_uniq_ambig = generate_salmon_uniq_ambig,
			gene_feature_gff_to_create_transcriptome_host = gene_feature_gff_to_create_transcriptome_host,
			gene_feature_gff_to_create_transcriptome_pathogen = gene_feature_gff_to_create_transcriptome_pathogen,
			gene_attribute_gff_to_create_transcriptome_host = gene_attribute_gff_to_create_transcriptome_host,
			gene_attribute_gff_to_create_transcriptome_pathogen = gene_attribute_gff_to_create_transcriptome_pathogen,
			run_salmon_selective_alignment = run_salmon_selective_alignment,
			kmer_length = kmer_length,
			writeUnmappedNames = writeUnmappedNames,
			softclipOverhangs = softclipOverhangs,
			dumpEq = dumpEq,
			writeMappings = writeMappings,
			keepDuplicates = keepDuplicates,
			salmon_sa_params_index = salmon_sa_params_index,
			salmon_sa_params_mapping = salmon_sa_params_mapping,
			run_salmon_alignment_based_mode = run_salmon_alignment_based_mode,
			salmon_alignment_based_params = salmon_alignment_based_params,
			run_star = run_star,
			outSAMunmapped = outSAMunmapped,
			outSAMattributes = outSAMattributes,
			outFilterMultimapNmax = outFilterMultimapNmax,
			outFilterType = outFilterType,
			alignSJoverhangMin = alignSJoverhangMin,
			alignSJDBoverhangMin = alignSJDBoverhangMin,
			outFilterMismatchNmax = outFilterMismatchNmax,
			outFilterMismatchNoverReadLmax = outFilterMismatchNoverReadLmax,
			alignIntronMin = alignIntronMin,
			alignIntronMax = alignIntronMax,
			alignMatesGapMax = alignMatesGapMax,
			limitBAMsortRAM = limitBAMsortRAM,
			winAnchorMultimapNmax = winAnchorMultimapNmax,
			sjdbOverhang = sjdbOverhang,
			quantTranscriptomeBan = quantTranscriptomeBan,
			star_salmon_index_params = star_salmon_index_params,
			star_salmon_alignment_params = star_salmon_alignment_params,
			outWigType = outWigType,
			outWigStrand = outWigStrand,
			star_index_params = star_index_params,
			star_alignment_params = star_alignment_params,
			run_htseq_uniquely_mapped = run_htseq_uniquely_mapped,
			stranded = stranded,
			max_reads_in_buffer = max_reads_in_buffer,
			minaqual = minaqual,
			htseq_params = htseq_params,
			gene_feature_gff_to_quantify_host = gene_feature_gff_to_quantify_host,
			host_gff_attribute = host_gff_attribute,
			gene_feature_gff_to_quantify_pathogen = gene_feature_gff_to_quantify_pathogen,
			pathogen_gff_attribute = pathogen_gff_attribute,
			mapping_statistics = mapping_statistics,
			rna_classes_to_replace_host = rna_classes_to_replace_host,
			publish_dir_mode = publish_dir_mode,
			email_on_fail = email_on_fail,
			plaintext_email = plaintext_email,
			max_multiqc_email_size = max_multiqc_email_size,
			monochrome_logs = monochrome_logs,
			multiqc_config = multiqc_config,
			tracedir = tracedir,
			help = help,
			email = email,
			custom_config_version = custom_config_version,
			custom_config_base = custom_config_base,
			hostnames = hostnames,
			config_profile_description = config_profile_description,
			config_profile_contact = config_profile_contact,
			config_profile_url = config_profile_url,
			outputbucket = thuuid.touchedbucket
            }
		output {
			Array[File] results = nfcoretask.results
		}
	}
task make_uuid {
	meta {
		volatile: true
}

command <<<
        python <<CODE
        import uuid
        print("gs://truwl-internal-inputs/nf-dualrnaseq/{}".format(str(uuid.uuid4())))
        CODE
>>>

  output {
    String uuid = read_string(stdout())
  }
  
  runtime {
    docker: "python:3.8.12-buster"
  }
}

task touch_uuid {
    input {
        String outputbucket
    }

    command <<<
        echo "sentinel" > sentinelfile
        gsutil cp sentinelfile ~{outputbucket}/sentinelfile
    >>>

    output {
        String touchedbucket = outputbucket
    }

    runtime {
        docker: "google/cloud-sdk:latest"
    }
}

task fetch_results {
    input {
        String outputbucket
        File execution_trace
    }
    command <<<
        cat ~{execution_trace}
        echo ~{outputbucket}
        mkdir -p ./resultsdir
        gsutil cp -R ~{outputbucket} ./resultsdir
    >>>
    output {
        Array[File] results = glob("resultsdir/*")
    }
    runtime {
        docker: "google/cloud-sdk:latest"
    }
}

task run_nfcoretask {
    input {
        String outputbucket
		String? name
		File samplesheet = "data/*{1,2}.fastq.gz"
		Boolean? single_end
		String outdir = "./results"
		Int max_cpus = 16
		String max_memory = "128.GB"
		String max_time = "240.h"
		String? fasta_host
		String? fasta_pathogen
		String? gff_host_genome
		String? gff_host_tRNA
		String? gff_pathogen
		String? transcriptome_host
		String? transcriptome_pathogen
		Boolean? read_transcriptome_fasta_host_from_file
		Boolean? read_transcriptome_fasta_pathogen_from_file
		String genome_host = "GRCh38"
		String genome_pathogen = "SL1344"
		Boolean? genomes_ignore
		Boolean? skip_fastqc
		String? fastqc_params
		Boolean? run_cutadapt
		String a = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
		String A = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
		Int quality_cutoff = 10
		String? cutadapt_params
		Boolean? run_bbduk
		Int minlen = 18
		String qtrim = "r"
		Int trimq = 10
		String ktrim = "r"
		Int k = 17
		Int mink = 11
		Int hdist = 1
		String adapters = "data/adapters.fa"
		String? bbduk_params
		String? libtype
		Int? incompatPrior
		Boolean? generate_salmon_uniq_ambig
		String gene_feature_gff_to_create_transcriptome_host = "['exon', 'tRNA']"
		String gene_feature_gff_to_create_transcriptome_pathogen = "['gene', 'sRNA', 'tRNA', 'rRNA']"
		String gene_attribute_gff_to_create_transcriptome_host = "transcript_id"
		String gene_attribute_gff_to_create_transcriptome_pathogen = "locus_tag"
		Boolean? run_salmon_selective_alignment
		Int kmer_length = 21
		Boolean? writeUnmappedNames
		Boolean? softclipOverhangs
		Boolean? dumpEq
		Boolean? writeMappings
		Boolean? keepDuplicates
		String? salmon_sa_params_index
		String? salmon_sa_params_mapping
		Boolean? run_salmon_alignment_based_mode
		String? salmon_alignment_based_params
		Boolean? run_star
		String outSAMunmapped = "Within"
		String outSAMattributes = "Standard"
		Int outFilterMultimapNmax = 999
		String outFilterType = "BySJout"
		Int alignSJoverhangMin = 8
		Int alignSJDBoverhangMin = 1
		Int outFilterMismatchNmax = 999
		Int outFilterMismatchNoverReadLmax = 1
		Int alignIntronMin = 20
		Int alignIntronMax = 1000000
		Int alignMatesGapMax = 1000000
		Int limitBAMsortRAM = 0
		Int winAnchorMultimapNmax = 999
		Int sjdbOverhang = 100
		String quantTranscriptomeBan = "Singleend"
		String? star_salmon_index_params
		String? star_salmon_alignment_params
		String outWigType = "None"
		String outWigStrand = "Stranded"
		String? star_index_params
		String? star_alignment_params
		Boolean? run_htseq_uniquely_mapped
		String stranded = "yes"
		Int max_reads_in_buffer = 30000000
		Int minaqual = 10
		String? htseq_params
		String gene_feature_gff_to_quantify_host = "['exon', 'tRNA']"
		String host_gff_attribute = "gene_id"
		String gene_feature_gff_to_quantify_pathogen = "['gene', 'sRNA', 'tRNA', 'rRNA']"
		String pathogen_gff_attribute = "locus_tag"
		Boolean? mapping_statistics
		String rna_classes_to_replace_host = "{base_dir}/data/RNA_classes_to_replace.csv"
		String publish_dir_mode = "copy"
		String? email_on_fail
		Boolean? plaintext_email
		String max_multiqc_email_size = "25.MB"
		Boolean? monochrome_logs
		String? multiqc_config
		String tracedir = "./results/pipeline_info"
		Boolean? help
		String? email
		String custom_config_version = "master"
		String custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/master"
		String? hostnames
		String? config_profile_description
		String? config_profile_contact
		String? config_profile_url

	}
	command <<<
		export NXF_VER=21.10.5
		export NXF_MODE=google
		echo ~{outputbucket}
		/nextflow -c /truwl.nf.config run /dualrnaseq-1.0.0  -profile truwl  --input ~{samplesheet} 	~{"--name " + name}	~{"--samplesheet " + samplesheet}	~{true="--single_end  " false="" single_end}	~{"--outdir " + outdir}	~{"--max_cpus " + max_cpus}	~{"--max_memory " + max_memory}	~{"--max_time " + max_time}	~{"--fasta_host " + fasta_host}	~{"--fasta_pathogen " + fasta_pathogen}	~{"--gff_host_genome " + gff_host_genome}	~{"--gff_host_tRNA " + gff_host_tRNA}	~{"--gff_pathogen " + gff_pathogen}	~{"--transcriptome_host " + transcriptome_host}	~{"--transcriptome_pathogen " + transcriptome_pathogen}	~{true="--read_transcriptome_fasta_host_from_file  " false="" read_transcriptome_fasta_host_from_file}	~{true="--read_transcriptome_fasta_pathogen_from_file  " false="" read_transcriptome_fasta_pathogen_from_file}	~{"--genome_host " + genome_host}	~{"--genome_pathogen " + genome_pathogen}	~{true="--genomes_ignore  " false="" genomes_ignore}	~{true="--skip_fastqc  " false="" skip_fastqc}	~{"--fastqc_params " + fastqc_params}	~{true="--run_cutadapt  " false="" run_cutadapt}	~{"--a " + a}	~{"--A " + A}	~{"--quality_cutoff " + quality_cutoff}	~{"--cutadapt_params " + cutadapt_params}	~{true="--run_bbduk  " false="" run_bbduk}	~{"--minlen " + minlen}	~{"--qtrim " + qtrim}	~{"--trimq " + trimq}	~{"--ktrim " + ktrim}	~{"--k " + k}	~{"--mink " + mink}	~{"--hdist " + hdist}	~{"--adapters " + adapters}	~{"--bbduk_params " + bbduk_params}	~{"--libtype " + libtype}	~{"--incompatPrior " + incompatPrior}	~{true="--generate_salmon_uniq_ambig  " false="" generate_salmon_uniq_ambig}	~{"--gene_feature_gff_to_create_transcriptome_host " + gene_feature_gff_to_create_transcriptome_host}	~{"--gene_feature_gff_to_create_transcriptome_pathogen " + gene_feature_gff_to_create_transcriptome_pathogen}	~{"--gene_attribute_gff_to_create_transcriptome_host " + gene_attribute_gff_to_create_transcriptome_host}	~{"--gene_attribute_gff_to_create_transcriptome_pathogen " + gene_attribute_gff_to_create_transcriptome_pathogen}	~{true="--run_salmon_selective_alignment  " false="" run_salmon_selective_alignment}	~{"--kmer_length " + kmer_length}	~{true="--writeUnmappedNames  " false="" writeUnmappedNames}	~{true="--softclipOverhangs  " false="" softclipOverhangs}	~{true="--dumpEq  " false="" dumpEq}	~{true="--writeMappings  " false="" writeMappings}	~{true="--keepDuplicates  " false="" keepDuplicates}	~{"--salmon_sa_params_index " + salmon_sa_params_index}	~{"--salmon_sa_params_mapping " + salmon_sa_params_mapping}	~{true="--run_salmon_alignment_based_mode  " false="" run_salmon_alignment_based_mode}	~{"--salmon_alignment_based_params " + salmon_alignment_based_params}	~{true="--run_star  " false="" run_star}	~{"--outSAMunmapped " + outSAMunmapped}	~{"--outSAMattributes " + outSAMattributes}	~{"--outFilterMultimapNmax " + outFilterMultimapNmax}	~{"--outFilterType " + outFilterType}	~{"--alignSJoverhangMin " + alignSJoverhangMin}	~{"--alignSJDBoverhangMin " + alignSJDBoverhangMin}	~{"--outFilterMismatchNmax " + outFilterMismatchNmax}	~{"--outFilterMismatchNoverReadLmax " + outFilterMismatchNoverReadLmax}	~{"--alignIntronMin " + alignIntronMin}	~{"--alignIntronMax " + alignIntronMax}	~{"--alignMatesGapMax " + alignMatesGapMax}	~{"--limitBAMsortRAM " + limitBAMsortRAM}	~{"--winAnchorMultimapNmax " + winAnchorMultimapNmax}	~{"--sjdbOverhang " + sjdbOverhang}	~{"--quantTranscriptomeBan " + quantTranscriptomeBan}	~{"--star_salmon_index_params " + star_salmon_index_params}	~{"--star_salmon_alignment_params " + star_salmon_alignment_params}	~{"--outWigType " + outWigType}	~{"--outWigStrand " + outWigStrand}	~{"--star_index_params " + star_index_params}	~{"--star_alignment_params " + star_alignment_params}	~{true="--run_htseq_uniquely_mapped  " false="" run_htseq_uniquely_mapped}	~{"--stranded " + stranded}	~{"--max_reads_in_buffer " + max_reads_in_buffer}	~{"--minaqual " + minaqual}	~{"--htseq_params " + htseq_params}	~{"--gene_feature_gff_to_quantify_host " + gene_feature_gff_to_quantify_host}	~{"--host_gff_attribute " + host_gff_attribute}	~{"--gene_feature_gff_to_quantify_pathogen " + gene_feature_gff_to_quantify_pathogen}	~{"--pathogen_gff_attribute " + pathogen_gff_attribute}	~{true="--mapping_statistics  " false="" mapping_statistics}	~{"--rna_classes_to_replace_host " + rna_classes_to_replace_host}	~{"--publish_dir_mode " + publish_dir_mode}	~{"--email_on_fail " + email_on_fail}	~{true="--plaintext_email  " false="" plaintext_email}	~{"--max_multiqc_email_size " + max_multiqc_email_size}	~{true="--monochrome_logs  " false="" monochrome_logs}	~{"--multiqc_config " + multiqc_config}	~{"--tracedir " + tracedir}	~{true="--help  " false="" help}	~{"--email " + email}	~{"--custom_config_version " + custom_config_version}	~{"--custom_config_base " + custom_config_base}	~{"--hostnames " + hostnames}	~{"--config_profile_description " + config_profile_description}	~{"--config_profile_contact " + config_profile_contact}	~{"--config_profile_url " + config_profile_url}	-w ~{outputbucket}
	>>>
        
    output {
        File execution_trace = "pipeline_execution_trace.txt"
        Array[File] results = glob("results/*/*html")
    }
    runtime {
        docker: "truwl/nfcore-dualrnaseq:1.0.0_0.1.0"
        memory: "2 GB"
        cpu: 1
    }
}
    