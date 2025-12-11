version 1.0

import "../wdl-common/wdl/structs.wdl"
import "../wdl-common/wdl/tasks/glnexus.wdl" as Glnexus
import "../wdl-common/wdl/tasks/sawfish.wdl" as Sawfish
import "../wdl-common/wdl/tasks/bcftools.wdl" as Bcftools

workflow joint {
  meta {
    description: "Tasks for joint-calling variants from a set of samples and splitting the joint calls by sample for parallel phasing."
  }

  parameter_meta {
    family_id: {
      name: "Cohort ID"
    }
    sample_ids: {
      name: "Sample IDs"
    }
    gvcfs: {
      name: "GVCFs"
    }   
    gvcf_indices: {
      name: "GVCF Indices"
    }

    aligned_bams: {
      name: "Aligned BAMs"
    }
    aligned_bam_indices: {
      name: "Aligned BAM Indices"
    }
    ref_map_file: {
      name: "Reference Map File"
    }
    glnexus_mem_gb: {
      name: "GLnexus Memory (GB)"
    }
    default_runtime_attributes: {
      name: "Default Runtime Attribute Struct"
    }

    split_joint_small_variant_vcfs: {
      name: "Joint-call small variant VCF, split by sample"
    }
    split_joint_small_variant_vcf_indices: {
      name: "Joint-call small variant VCF indices, split by sample"
    }
    sv_supporting_reads: {
      name: "Supporting reads JSON"
    }
    sv_copynum_bedgraph: {
      name: "Copy number bedgraph"
    }
    sv_depth_bw: {
      name: "Depth bigWig"
    }
    sv_gc_bias_corrected_depth_bw: {
      name: "GC bias corrected depth bigWig"
    }
    sv_maf_bw: {
      name: "MAF bigWig"
    }
    sv_copynum_summary: {
      name: "Copy number summary JSON"
    }
  }

  input {
    String family_id
    Array[String] sample_ids

    Array[File] gvcfs
    Array[File] gvcf_indices


    Array[File] aligned_bams
    Array[File] aligned_bam_indices

    File ref_map_file

    Int? glnexus_mem_gb

    RuntimeAttributes default_runtime_attributes
  }

  Map[String, String] ref_map = read_map(ref_map_file)



  call Glnexus.glnexus {
    input:
      cohort_id          = family_id + ".joint",
      gvcfs              = gvcfs,
      gvcf_indices       = gvcf_indices,
      ref_name           = ref_map["name"],
      mem_gb             = glnexus_mem_gb,
      runtime_attributes = default_runtime_attributes
  }

  String glnexus_vcf_basename = basename(glnexus.vcf, ".vcf.gz")

  scatter (sample_id in sample_ids) {
    String split_glnexus_vcf_name = "~{sample_id}.~{glnexus_vcf_basename}.vcf.gz"
    String split_glnexus_vcf_index_name = "~{sample_id}.~{glnexus_vcf_basename}.vcf.gz.tbi"
  
    call Bcftools.split_vcf_by_sample as split_glnexus {
      input:
        sample_id             = sample_id,
        vcf                   = glnexus.vcf,
        vcf_index             = glnexus.vcf_index,
        split_vcf_name        = split_glnexus_vcf_name,
        split_vcf_index_name  = split_glnexus_vcf_index_name,
        runtime_attributes    = default_runtime_attributes
    }
  }

  output {
    Array[File] split_joint_small_variant_vcfs             = split_glnexus.split_vcf
    Array[File] split_joint_small_variant_vcf_indices      = split_glnexus.split_vcf_index
  }
}
