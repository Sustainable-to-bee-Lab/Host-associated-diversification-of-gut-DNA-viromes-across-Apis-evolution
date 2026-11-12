# Apis gut DNA virome core analysis code

## Contents

- [01 Viral identification, validation, and taxonomic assignment](#module-01)
- [02 Definition of viral units and viral clustering](#module-02)
- [03 Diversity analyses and phage-bacterium associations](#module-03)
- [04 Bacterial host prediction for bacteriophages](#module-04)
- [05 Prediction and functional annotation of viral proteins](#module-05)
- [06 Phylogenetic analysis](#module-06)
- [07 Phage-bacterium network construction and robustness analysis](#module-07)


<a id="module-01"></a>

## 01 Viral identification, validation, and taxonomic assignment

<a id="detection-and-filtering-with-five-methods"></a>

### Detection and filtering with five methods

<!-- CODE detection-and-filtering-with-five-methods BEGIN -->
````bash
while IFS= read -r sample; do
    assembly="${sample}.len500.fna"
    diamond blastx --threads 20 -d "$NR_DB" -q "$assembly" \
        -f 6 staxids qseqid sseqid pident length qlen qstart qend slen sstart send evalue bitscore \
        -o "${sample}.nr_blastx.tsv" --sensitive -k 1 -e 1e-5 --log

    genomad end-to-end --cleanup "$assembly" "${sample}.genomad" "$GENOMAD_DB"
    awk -F '\t' 'NR>1 && $7>0.8 {print $1}' "$GENOMAD_SUMMARY" > "${sample}.genomad.id"

    virsorter run --keep-original-seq -i "$assembly" -w "${sample}.virsorter2" \
        --min-length 5000 --min-score 0.7 -j 28 all

    viralverify -f "$assembly" -o "${sample}.viralverify" --hmm "$VIRALVERIFY_HMM"
    seqkit seq -m 5001 "$VIRALVERIFY_VIRUS_FASTA" > "${sample}.viralverify.gt5000.fa"
    checkv end_to_end "${sample}.viralverify.gt5000.fa" "${sample}.viralverify.checkv" \
        -t 28 -d "$CHECKV_DB"

    python "$DVF_SCRIPT" -i "$assembly" -o "${sample}.dvf" -l 5000 -c 20
    awk -F '\t' 'NR>1 && $2>5000 && $3>0.9 && $4<0.01 {print $1}' \
        "$DVF_PREDICTIONS" > "${sample}.dvf.id"
done < samples.txt
````
<!-- CODE detection-and-filtering-with-five-methods END -->

<a id="candidate-quality-assessment-and-nt-alignment"></a>

### Candidate quality assessment and nt alignment

<!-- CODE candidate-quality-assessment-and-nt-alignment BEGIN -->
````bash
checkv end_to_end candidates.fa checkv_candidates -t 30 -d "$CHECKV_DB"
blastn -db "$NT_DB" -query candidates.fa \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids' \
    -out candidates.blastn.tsv -evalue 1e-5
taxonkit lineage --data-dir "$TAXONOMY_DB" taxids.txt > taxid_lineage.tsv
````
<!-- CODE candidate-quality-assessment-and-nt-alignment END -->

<a id="viral-confidence-classification"></a>

### Viral-confidence classification

<!-- CODE viral-confidence-classification BEGIN -->
````bash
awk 'BEGIN { FS=OFS="\t" }
{
    num_methods=split($14, methods, " ")
    if (num_methods>1 || $12=="Complete" || $12=="High-quality" ||
        $12=="Medium-quality" || $5=="Virus") {
        confidence="High"
    } else if ($12=="Low-quality" || $12=="Not-determined") {
        confidence=($3=="unknown" || $3+0<60) ? "Medium" : "Low"
    } else confidence="Low"
    print $0,confidence
}' evidence.tsv > confidence.tsv
````
<!-- CODE viral-confidence-classification END -->

<a id="protein-evidence-filtering"></a>

### Protein-evidence filtering

<!-- CODE protein-evidence-filtering BEGIN -->
````bash
awk -F '\t' 'BEGIN {OFS="\t"}
NR==1 {print $0,"BlastxResult";next}
{
    confidence=$15;identity=$16+0
    protein=tolower($19);tax=tolower($20);result=""
    if (confidence=="High" || confidence=="Medium") {
        if (protein ~ /(virus|phage|viria|virea|viricota|viricetes|virales|viridae|rna-dependent rna polymerase|rdrp)/ ||
            tax ~ /(virus|phage|viria|virea|viricota|viricetes|virales|viridae|rna-dependent rna polymerase|rdrp)/)
            result="virus protein"
        else if (protein ~ /(hypothetical|uncharacterized)/ || tax ~ /(hypothetical|uncharacterized)/)
            result="hypothetical/uncharacterized Protein"
        else if (identity<60 || protein=="" || tax=="") result="unknown protein"
        else result="others protein"
    }
    print $0,result
}' confidence_blastx.tsv > classified.tsv

awk -F '\t' '$21=="virus protein" || $21=="hypothetical/uncharacterized Protein" ||
    $21=="unknown protein" {print $1}' classified.tsv > final_virus.id
seqkit grep -f final_virus.id candidates.fa > retained_virus.fa
````
<!-- CODE protein-evidence-filtering END -->

<a id="prophage-integration-and-final-length-filtering"></a>

### Prophage integration and final length filtering

BED regions are extracted from geNomad and VirSorter2 and use zero-based, half-open coordinates.

<!-- CODE prophage-integration-and-final-length-filtering BEGIN -->
````bash
sort -k1,1 -k2,2n prophage_regions.bed | bedtools merge -i - > prophage_regions.merged.bed
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1"|provirus_"($2+1)"_"$3}' \
    prophage_regions.merged.bed > prophage_regions.named.bed
bedtools getfasta -fi parent_contigs.fa -bed prophage_regions.named.bed \
    -name -fo integrated_prophage.raw.fa
seqkit replace -p '::.*$' -r '' integrated_prophage.raw.fa > integrated_prophage.fa
cut -f1 prophage_regions.merged.bed | sort -u > prophage_parent.id
seqkit grep -v -f prophage_parent.id retained_virus.fa > retained_without_parent.fa
cat integrated_prophage.fa retained_without_parent.fa > merged_virus.fa
seqkit rmdup -s merged_virus.fa | seqkit seq -m 1000 > final_virus.min1000.fa
````
<!-- CODE prophage-integration-and-final-length-filtering END -->

<a id="soap-indexing-mapping-and-coverage"></a>

<a id="final-taxonomic-integration"></a>

### Final taxonomic integration

<!-- CODE final-taxonomic-integration BEGIN -->
````r
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("Usage: Rscript script.R input_bundle.rds output.rds")
bundle <- readRDS(args[1])
required <- c("meta","phabox_lineage","previous_taxonomy")
stopifnot(is.list(bundle), all(required %in% names(bundle)))
list2env(bundle[required], envir=environment())
stopifnot(identical(names(meta), c("sample","Host","location")),
  nrow(meta)==69L, !anyNA(meta), !anyDuplicated(meta$sample),
  setequal(meta$location,c("JL_DH","YN_HH","YN_JH","YN_NJ","YN_PE","YN_YX")))
informative_taxonomy <- function(x) {
  x <- x[!is.na(x) & !x %in% c("","-","Viruses","Unclassified","no hits to database")]
  if(length(x))x[1]else NA_character_
}
taxonomy <- phabox_lineage |>
  dplyr::rowwise() |>
  dplyr::mutate(new_taxonomy=informative_taxonomy(c(family,order,class))) |>
  dplyr::ungroup() |>
  dplyr::left_join(previous_taxonomy,by="contig_id") |>
  dplyr::mutate(taxonomy=dplyr::coalesce(new_taxonomy,old_taxonomy,"Unclassified"))
stopifnot(!anyDuplicated(taxonomy$contig_id))
saveRDS(taxonomy,args[2])
````
<!-- CODE final-taxonomic-integration END -->


### Split final sequences by sample, build SOAP indices, and map reads

<!-- CODE soap-indexing-mapping-and-coverage BEGIN -->
````bash
seqkit replace -p ' .*$' -r '' final_virus.min1000.fa > final_virus.clean.fa
while IFS= read -r sample; do
    seqkit grep -r -p "^${sample}_k[0-9]+_" final_virus.clean.fa > "${sample}.final_virus.fa"
    2bwt-builder "${sample}.final_virus.fa"
    soap -a "${sample}.clean.1.fq.gz" -b "${sample}.clean.2.fq.gz" \
        -D "${sample}.final_virus.fa.index" \
        -o "${sample}.assembly.aln" -2 "${sample}.assembly.se.aln" \
        -M 4 -l 30 -r 1 -v 3 -m 200 -p 12
    soap.coverage -cvg -i "${sample}.assembly.aln" "${sample}.assembly.se.aln" \
        -refsingle "${sample}.final_virus.fa" -o "${sample}.result.txt"
done < samples.txt
````
<!-- CODE soap-indexing-mapping-and-coverage END -->

<a id="extract-cvg-and-depth-from-soap-results"></a>

### Extract depth from SOAP results

<!-- CODE extract-cvg-and-depth-from-soap-results BEGIN -->
````bash
while IFS= read -r sample; do
    perl - "${sample}.result.txt" "$sample" <<'PERL'
use strict;
use warnings;
my ($result_file,$sample)=@ARGV;
open my $in,'<',$result_file or die "$result_file: $!";
open my $dep,'>',"$sample.depth" or die $!;
while (<$in>) {
    chomp;
    next if $_ eq '';
    my ($id)=split /:/;
    my @a=split /\s+|\//;
    my $depth=$a[-1];
    $depth =~ s/Depth://;
    next if $depth eq '' || $depth !~ /^[0-9.eE+-]+$/;
    print {$dep} "$id\t$depth\n";
}
PERL
done < samples.txt
````
<!-- CODE extract-cvg-and-depth-from-soap-results END -->


<a id="module-02"></a>

## 02 Definition of viral units and viral clustering

<a id="vcontact2-clustering"></a>

### vConTACT2 clustering

<!-- CODE vcontact2-clustering BEGIN -->
````bash
for branch in phage other; do
    vcontact2_gene2genome --proteins "${branch}.faa" --output "${branch}.gene2genome.csv" -s Prodigal-FAA
    vcontact2 --rel-mode Diamond --raw-proteins "${branch}.faa" \
        --proteins-fp "${branch}.gene2genome.csv" --db ProkaryoticViralRefSeq207-Merged \
        --output-dir "${branch}.vcontact2_result"
done
````
<!-- CODE vcontact2-clustering END -->

<a id="host-pair-vc-sharing-and-exact-permutation"></a>


<a id="rebuild-branch-units-and-depth-abundance-from-current-vcs"></a>

### Rebuild branch units and Depth abundance from current VCs

<!-- CODE rebuild-branch-units-and-depth-abundance-from-current-vcs BEGIN -->
````python
import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


def norm(seq):
    return seq.replace('|provirus_', '_provirus_')


def read_table(path, delimiter=','):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_rows(path, rows, delimiter='\t'):
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        csv.writer(handle, delimiter=delimiter, lineterminator='\n').writerows(rows)


def run_perl(script, *args):
    subprocess.run(['perl', str(script)] + [str(x) for x in args],
                   check=True, stdout=subprocess.DEVNULL)


def main():
    parser = argparse.ArgumentParser()
    for name in ('samples', 'catalog', 'phage-taxonomy', 'other-taxonomy',
                 'phage-vcs', 'other-vcs', 'depth-dir', 'coverage-dir', 'abundance-scripts', 'out'):
        parser.add_argument('--' + name, required=True, type=Path)
    args = parser.parse_args()
    samples = args.samples.read_text(encoding='utf-8-sig').split()
    catalog = read_table(args.catalog)
    final = {norm(row['contig_id']): row for row in catalog}
    assert len(samples) == len(set(samples)) == 69
    assert len(catalog) == len(final) == 13661
    assert {r['Sample'] for r in catalog} == set(samples)
    assert all(int(r['contig_length']) >= 1000 for r in catalog)
    args.out.mkdir(parents=True, exist_ok=False)
    work = args.out / 'contig_depth'
    work.mkdir()
    source_paths = [args.samples, args.catalog, args.phage_taxonomy,
                    args.other_taxonomy, args.phage_vcs, args.other_vcs]
    allalpha = args.abundance_scripts / 'AllAlpha.pl'
    ko_stat = args.abundance_scripts / 'KO_stat.pl'
    pct = args.abundance_scripts / 'pct.pl'
    source_paths.extend([allalpha, ko_stat, pct])
    depth_ids = {}
    zero_coverage_rows = []
    for sample in samples:
        depth_file = args.depth_dir / (sample + '.depth')
        source_paths.append(depth_file)
        ids = {}
        with depth_file.open() as handle:
            for line in handle:
                fields = line.split()
                if not fields:
                    continue
                seq, value = fields[:2]
                if seq == 'sample':
                    continue
                key = norm(seq)
                assert key not in ids
                assert math.isfinite(float(value)) and float(value) >= 0
                ids[key] = seq
        needed = {key for key, row in final.items() if row['Sample'] == sample}
        missing = needed - set(ids)
        depth_input = depth_file
        if missing:
            coverage_file = args.coverage_dir / (sample + '.result.txt')
            source_paths.append(coverage_file)
            recovered = []
            with coverage_file.open() as handle:
                for line in handle:
                    record = re.match(r'^(.+?):\s*(\d+)/(\d+)\s+.*Depth:([^\s]+)', line)
                    if not record or norm(record.group(1)) not in missing:
                        continue
                    seq, covered, length, depth = record.groups()
                    key = norm(seq)
                    assert int(covered) == 0 and depth.lower() in ('nan', '0', '0.0')
                    assert int(length) == int(final[key]['contig_length']) and key not in ids
                    ids[key] = seq
                    recovered.append([seq, '0'])
                    zero_coverage_rows.append([sample, key, int(length), 0, depth])
            assert needed <= set(ids), (sample, 'missing final depth records')
            depth_input = work / (sample + '.complete.depth')
            with depth_file.open() as handle:
                existing = list(csv.reader(handle, delimiter='\t'))
            write_rows(depth_input, existing + recovered)
        depth_ids[sample] = ids
        run_perl(allalpha, work / sample, depth_input)
    write_rows(args.out / 'zero_coverage_depth_records.tsv',
               [['Sample', 'sequence_ID', 'Length_bp', 'Covered_bp', 'SOAP_Depth']] + zero_coverage_rows)

    summary, mapping = [], []
    all_selected = set()
    for branch, tax_path, vc_path, stem, expected_vcs in [
        ('phage', args.phage_taxonomy, args.phage_vcs, 'Apis_Phage', 1100),
        ('other', args.other_taxonomy, args.other_vcs, 'Apis_Virus', 773),
    ]:
        taxonomy = {}
        with tax_path.open() as handle:
            for fields in csv.reader(handle, delimiter='\t'):
                if not fields:
                    continue
                key = norm(fields[0])
                assert key not in taxonomy and key in final
                value = fields[1].strip()
                taxonomy[key] = 'Unclassified' if value in ('', 'NA', 'NaN', 'None') else value
                assert taxonomy[key] == final[key]['family_classification']
        assert not (all_selected & set(taxonomy))
        all_selected.update(taxonomy)
        network = {}
        for row in read_table(vc_path):
            key = norm(row['Genome'])
            assert key not in network
            network[key] = row
        groups = defaultdict(list)
        counts = defaultdict(Counter)
        for seq, tax in taxonomy.items():
            vc = network.get(seq, {}).get('VC', '')
            if vc:
                groups[vc].append(seq)
                counts[vc][tax] += 1
        assert len(groups) == expected_vcs
        consensus, rule_rows = {}, []
        for vc in sorted(groups):
            tax, count = sorted(counts[vc].items(), key=lambda x: (-x[1], x[0]))[0]
            total = len(groups[vc])
            use = count * 5 >= total * 4
            if use:
                consensus[vc] = re.sub(r'\s+', '_', tax) + '_' + vc
            rule_rows.append([vc, tax, count, total, count / total,
                              'majority_ge_80' if use else 'no_majority'])
        folder = args.out / ('phage_VCrelabun' if branch == 'phage' else 'virus_VCrelabun')
        index_dir = folder / 'family_VC_index_by_sample'
        profile_dir = folder / 'sample_profiles'
        index_dir.mkdir(parents=True)
        profile_dir.mkdir()
        per_sample = {s: [] for s in samples}
        branch_rows = []
        for seq, tax in taxonomy.items():
            sample = final[seq]['Sample']
            raw = network.get(seq, {})
            vc = raw.get('VC', '')
            unit = consensus.get(vc, tax)
            assert not re.search(r'\s', unit)
            suffix = re.search(r'(VC_\d+_\d+)$', unit)
            assert not suffix or suffix.group(1) == vc
            per_sample[sample].append([depth_ids[sample][seq], unit])
            row = [seq, sample, branch, tax, vc, raw.get('VC Status', ''), unit,
                   'majority_ge_80' if vc in consensus else 'original_taxonomy']
            branch_rows.append(row)
            mapping.append(row)
        headers = ['sequence_ID', 'Sample', 'Analysis_branch', 'Taxonomy_label',
                   'VC_id', 'VC_status', 'Taxonomy_VC_unit', 'Unit_rule']
        write_rows(folder / 'sequence_VC_unit_mapping.tsv', [headers] + branch_rows)
        write_rows(folder / (stem + '.family.index.renamed.tsv'),
                   [row for sample in samples for row in per_sample[sample]])
        write_rows(folder / 'vc_taxonomy_majority80.tsv',
                   [['VC', 'Major_taxonomy', 'Major_count', 'Total_count', 'Fraction', 'Rule']] + rule_rows)
        no_header_files = []
        for sample in samples:
            index = index_dir / (sample + '_family_VC.index')
            write_rows(index, per_sample[sample])
            prefix = profile_dir / sample
            run_perl(ko_stat, index, work / (sample + '.alpha.txt'), prefix)
            ko_file = profile_dir / (sample + '.KO.Abu.txt')
            cleaned = profile_dir / (sample + '.KO.Abu.noheader.txt')
            with ko_file.open() as handle:
                rows = list(csv.reader(handle, delimiter='\t'))
            write_rows(cleaned, [r for r in rows if r and r[0] != 'sample'])
            no_header_files.append(cleaned)
        prefix = folder / (stem + '_VC_depth')
        run_perl(allalpha, prefix, *no_header_files)
        alpha = folder / (stem + '_VC_depth.alpha.txt')
        with alpha.open() as handle:
            alpha_rows = list(csv.reader(handle, delimiter='\t'))
        assert alpha_rows[0][1:] == [str(x) for x in no_header_files]
        alpha_rows[0] = ['sample'] + samples
        assert len({r[0] for r in alpha_rows[1:]}) == len(alpha_rows) - 1
        assert {r[0] for r in alpha_rows[1:]} == {r[6] for r in branch_rows}
        write_rows(alpha, alpha_rows)
        percentage = folder / (stem + '_VC_depth.alpha.pct.txt')
        run_perl(pct, alpha, percentage)
        with percentage.open() as handle:
            pct_rows = list(csv.reader(handle, delimiter='\t'))
        assert pct_rows[0] == ['sample'] + samples
        sums = [sum(float(r[j]) for r in pct_rows[1:]) for j in range(1, 70)]
        assert max(abs(x - 100) for x in sums) < 1e-8
        write_rows(folder / (stem + '_VC_depth.alpha.pct.csv'), pct_rows, ',')
        summary.append({'branch': branch, 'sequences': len(taxonomy),
                        'raw_VCs': len(groups), 'consensus_VCs': len(consensus),
                        'features': len(pct_rows) - 1, 'samples': len(samples),
                        'index_files': len(per_sample),
                        'missing_network_rows': sum(s not in network for s in taxonomy),
                        'max_percentage_sum_error': max(abs(x - 100) for x in sums)})
    assert all_selected == set(final)
    write_rows(args.out / 'sequence_VC_unit_mapping.tsv', [headers] + mapping)
    manifest = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    result = {'samples': 69, 'sequences': 13661, 'branches': summary,
              'consensus_threshold': 0.8, 'consensus_denominator': 'final study contigs in each VC',
              'abundance': 'SOAP Depth summed within unit; AllAlpha one-decimal rounding at contig and unit stages; separate scope percentages',
              'source_sha256': manifest, 'zero_coverage_depth_records': len(zero_coverage_rows),
              'downstream_analyses_rerun': False}
    (args.out / 'summary.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
````
<!-- CODE rebuild-branch-units-and-depth-abundance-from-current-vcs END -->


### VC sharing: shared matrices, distances and counts

<!-- CODE host-pair-vc-sharing-and-exact-permutation BEGIN -->
````r
# 1. Necessary inputs and checks / 必要输入与检查
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else "."
out <- if(length(args)>1) args[2] else "vc_preparation"
samples <- scan(file.path(root,"samples.txt"),what="",quiet=TRUE)
meta <- read.csv(file.path(root,"metadata_69.csv"),colClasses="character")
info <- read.csv(file.path(root,"Apis_all_info/Apis_all_info/Apis_virome_all_info.csv"),
  colClasses="character",check.names=FALSE,na.strings=c("","NA"))
stopifnot(length(samples)==69L,!anyDuplicated(samples),!anyDuplicated(meta$sample),
  setequal(samples,meta$sample),nrow(info)==13661L,!anyDuplicated(info$contig_id))
meta <- meta[match(samples,meta$sample),]
hosts <- c("Aand","Aflo","Acer","Amel","Ador","Alab")
subgenus <- c(Aand="dwarf",Aflo="dwarf",Acer="cavity",Amel="cavity",Ador="giant",Alab="giant")
stopifnot(setequal(unique(meta$Host),hosts),identical(info$Sample,sub("_.*","",info$contig_id)),
  setequal(unique(info$Sample),samples),all(as.numeric(info$contig_length)>=1000),
  all(info$Host==meta$Host[match(info$Sample,meta$sample)]),
  !anyNA(meta[c("Host","location")]),setequal(meta$location,c("JL_DH","YN_HH","YN_JH","YN_NJ","YN_PE","YN_YX")))
valid_vc <- function(x) !is.na(x) & nzchar(x)
stopifnot(all(info$VC_scope[valid_vc(info$VC)] %in% c("phage","other")))
datasets <- lapply(c("phage","other"),function(g) info[which(info$VC_scope==g),])
names(datasets) <- c("corrected_phage","corrected_other")
bind <- function(x) { z <- do.call(rbind,x);rownames(z) <- NULL;z }

# 2. VC presence matrices and Jaccard distances / VC 矩阵和距离
presence_matrices <- shared_matrices <- union_matrices <- distance_matrices <- list()
for(label in names(datasets)) {
  x <- datasets[[label]];x <- x[valid_vc(x$VC),]
  units <- sort(unique(x$VC))
  presence <- matrix(0L,69L,length(units),dimnames=list(samples,units))
  presence[cbind(match(x$Sample,samples),match(x$VC,units))] <- 1L
  stopifnot(all(rowSums(presence)>0))
  shared <- tcrossprod(presence)
  union <- outer(rowSums(presence),rowSums(presence),"+")-shared
  presence_matrices[[label]] <- presence
  shared_matrices[[label]] <- shared
  union_matrices[[label]] <- union
  distance_matrices[[label]] <- 1-shared/union
}
sample_index <- t(combn(seq_along(samples),2))
host_code <- match(meta$Host,hosts)
subgenus_code <- match(subgenus[hosts],unique(subgenus[hosts]))
pair_category <- function(labels,index) {
  left <- labels[index[,1]];right <- labels[index[,2]]
  ifelse(left==right,1L,ifelse(subgenus_code[left]==subgenus_code[right],2L,3L))
}
category_names <- c("Same host","Same subgenus","Different subgenera")
sample_category <- pair_category(host_code,sample_index)
sample_type <- category_names[sample_category]
stopifnot(identical(as.integer(tabulate(sample_category,nbins=3)),c(377L,382L,1587L)))
summaries <- categories <- host_counts <- sample_pairs <- list()
sample_sets <- split(meta$sample,factor(meta$Host,levels=hosts))
for(label in names(datasets)) {
  x <- datasets[[label]]
  vcx <- x[valid_vc(x$VC),]
  vc_hosts <- lapply(split(vcx$Host,vcx$VC),function(h) sort(unique(h)))
  cls <- vapply(vc_hosts,function(h) if(length(h)==1L) "Host-specific" else
    if(length(unique(subgenus[h]))==1L) "Subgenus-shared" else "Cross-subgenus",character(1))
  categories[[label]] <- data.frame(dataset=label,VC=names(vc_hosts),
    hosts=vapply(vc_hosts,paste,collapse=";",character(1)),n_hosts=lengths(vc_hosts),class=cls)
  summaries[[label]] <- data.frame(dataset=label,n_contigs=nrow(x),
    clustered_contigs=nrow(vcx),unclustered_contigs=nrow(x)-nrow(vcx),
    n_VC=length(vc_hosts),host_specific=sum(cls=="Host-specific"),
    subgenus_shared=sum(cls=="Subgenus-shared"),cross_subgenus=sum(cls=="Cross-subgenus"),
    all_six=sum(lengths(vc_hosts)==6L))
  sets <- lapply(hosts,function(h) unique(vcx$VC[vcx$Host==h]))
  host_counts[[label]] <- data.frame(dataset=label,Host=hosts,n_samples=lengths(sample_sets),
    n_contigs=vapply(hosts,function(h) sum(x$Host==h),integer(1)),n_VC=lengths(sets),
    n_specific=vapply(hosts,function(h) sum(cls=="Host-specific" &
      vapply(vc_hosts,function(hh) h %in% hh,logical(1))),integer(1)))
  value <- distance_matrices[[label]][sample_index]
  sample_pairs[[label]] <- data.frame(dataset=label,
    sample1=samples[sample_index[,1]],sample2=samples[sample_index[,2]],
    host1=meta$Host[sample_index[,1]],host2=meta$Host[sample_index[,2]],
    pair_type=sample_type,shared=shared_matrices[[label]][sample_index],
    union=union_matrices[[label]][sample_index],Jaccard=1-value,Jaccard_distance=value)
}
dir.create(out,recursive=TRUE,showWarnings=FALSE)
results <- list(vc_category_summary.csv=bind(summaries),vc_host_membership.csv=bind(categories),
  host_richness_and_sample_counts.csv=bind(host_counts),sample_pairwise_Jaccard.csv=bind(sample_pairs))
for(name in names(results)) write.csv(results[[name]],file.path(out,name),row.names=FALSE)
saveRDS(list(samples=samples,presence=presence_matrices,distance=distance_matrices),
  file.path(out,"VC_matrices.rds"))
````
<!-- CODE host-pair-vc-sharing-and-exact-permutation END -->


<a id="within-location-vc-sharing-pair-comparisons"></a>

### Within-location comparisons of VC sharing within and between hosts

<!-- CODE within-location-vc-sharing-pair-comparisons BEGIN -->
````r
if(!requireNamespace("permute",quietly=TRUE)) stop("The permute package is required.")
args <- commandArgs(trailingOnly=TRUE)
pair_file <- if(length(args)) args[1] else "sample_pairwise_Jaccard.csv"
metadata_file <- if(length(args)>1) args[2] else "metadata_69.csv"
out <- if(length(args)>2) args[3] else "within_location_tests"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
pairs <- read.csv(pair_file,check.names=FALSE,stringsAsFactors=FALSE)
meta <- read.csv(metadata_file,check.names=FALSE,stringsAsFactors=FALSE)
locations <- c("JL_DH","YN_JH","YN_YX")
datasets <- c("corrected_phage","corrected_other")
stopifnot(nrow(meta)==69L,!anyDuplicated(meta$sample),nrow(pairs)==4692L,
          setequal(unique(pairs$dataset),datasets),
          setequal(meta$location,c("JL_DH","YN_HH","YN_JH","YN_NJ","YN_PE","YN_YX")),
          all(is.finite(pairs$Jaccard_distance)),
          all(pairs$Jaccard_distance>=0 & pairs$Jaccard_distance<=1),
          max(abs(pairs$Jaccard_distance+pairs$Jaccard-1))<1e-12)
# Paired distances remain fixed; entire sample Host labels are reassigned.
# permute generates within-location sample-index permutations.
# For small reference sets, retain unique Host allocations until coverage is complete.
# The two-sided definition matches the existing pair-category analysis:
# twice the smaller inclusive permutation tail, capped at one.
tests <- summaries <- host_counts <- list()
for(li in seq_along(locations)) {
  loc <- locations[li]
  m <- meta[meta$location==loc,,drop=FALSE]
  labels <- as.integer(factor(m$Host))
  n <- nrow(m)
  idx <- t(combn(seq_len(n),2))
  same <- labels[idx[,1]]==labels[idx[,2]]
  counts <- tabulate(labels)
  n_distinct <- round(exp(lgamma(n+1)-sum(lgamma(counts+1))))
  exact <- n_distinct<=9999
  seed <- 20261105L+li
  design <- permute::how(within=permute::Within(type="free"),
    blocks=factor(m$location),nperm=9999L,complete=FALSE,
    minperm=0L,make=FALSE)
  set.seed(seed)
  draw_indices <- function() as.matrix(permute::shuffleSet(n,nset=9999L,
    control=design,check=FALSE))
  index_batches <- list(draw_indices())
  map_labels <- function(index) matrix(labels[index],nrow=nrow(index),ncol=n)
  perm_labels <- map_labels(index_batches[[1]])
  if(exact) {
    perm_labels <- unique(perm_labels)
    while(nrow(perm_labels)<n_distinct) {
      if(length(index_batches)>=100L) stop("Incomplete exact Host-label coverage.")
      index_batches[[length(index_batches)+1L]] <- draw_indices()
      perm_labels <- unique(rbind(perm_labels,map_labels(index_batches[[length(index_batches)]])))
    }
    stopifnot(nrow(perm_labels)==n_distinct)
    perm_labels <- perm_labels[do.call(order,as.data.frame(perm_labels)),,drop=FALSE]
  }
  permutation_indices <- do.call(rbind,index_batches)
  stopifnot(all(apply(permutation_indices,1L,function(i) identical(sort(i),seq_len(n)))),
    all(matrix(m$location[permutation_indices],nrow=nrow(permutation_indices))==
        matrix(m$location,nrow=nrow(permutation_indices),ncol=n,byrow=TRUE)),
    all(vapply(seq_along(counts),function(g) all(rowSums(perm_labels==g)==counts[g]),logical(1))))
  same_perm <- perm_labels[,idx[,1],drop=FALSE]==perm_labels[,idx[,2],drop=FALSE]
  stopifnot(all(rowSums(same_perm)==sum(same)))
  values <- matrix(NA_real_,nrow(idx),length(datasets),dimnames=list(NULL,datasets))
  for(k in seq_along(datasets)) {
    d <- pairs[pairs$dataset==datasets[k] & pairs$sample1 %in% m$sample &
               pairs$sample2 %in% m$sample,,drop=FALSE]
    a <- match(d$sample1,m$sample); b <- match(d$sample2,m$sample)
    stopifnot(!anyNA(c(a,b)),nrow(d)==nrow(idx),all(a!=b))
    mat <- matrix(NA_real_,n,n)
    mat[cbind(a,b)] <- d$Jaccard_distance
    mat[cbind(b,a)] <- d$Jaccard_distance
    values[,k] <- mat[idx]
    stopifnot(all(is.finite(values[,k])))
    for(g in c(TRUE,FALSE)) {
      v <- values[same==g,k]
      summaries[[length(summaries)+1L]] <- data.frame(dataset=datasets[k],
        location=loc,pair_type=if(g) "Same host" else "Different hosts",
        n_pairs=length(v),mean_distance=mean(v),median_distance=median(v),
        q25=unname(quantile(v,.25)),q75=unname(quantile(v,.75)),
        mean_sharing=mean(1-v))
    }
  }
  observed <- colMeans(values[!same,,drop=FALSE])-colMeans(values[same,,drop=FALSE])
  sum_same <- same_perm %*% values
  sum_different <- matrix(colSums(values),nrow(same_perm),ncol(values),byrow=TRUE)-sum_same
  null <- sum_different/sum(!same)-sum_same/sum(same)
  for(k in seq_along(datasets)) {
    upper <- sum(null[,k]>=observed[k]-1e-12)
    lower <- sum(null[,k]<=observed[k]+1e-12)
    if(exact) {
      p_upper <- upper/nrow(null); p_lower <- lower/nrow(null)
    } else {
      p_upper <- (upper+1)/(nrow(null)+1); p_lower <- (lower+1)/(nrow(null)+1)
    }
    tests[[length(tests)+1L]] <- data.frame(dataset=datasets[k],location=loc,
      n_samples=n,n_host_species=length(counts),same_host_pairs=sum(same),
      different_host_pairs=sum(!same),
      mean_same_host_distance=mean(values[same,k]),
      mean_different_host_distance=mean(values[!same,k]),
      statistic=observed[k],
      p_two_sided=min(1,2*min(p_upper,p_lower)),
      permutation_method=if(exact) "exact_unique_Host_allocations" else "Monte_Carlo_Host_labels",
      null_allocations=nrow(null),possible_distinct_allocations=n_distinct,
      unique_allocations_evaluated=nrow(unique(perm_labels)),seed=if(exact) NA_integer_ else seed,
      unit="whole_sample_Host_label",alternative="two_sided",
      allocation_generator="permute::shuffleSet",
      permute_version=as.character(utils::packageVersion("permute")),
      generation_seed=seed,index_permutations_generated=nrow(permutation_indices))
  }
  host_counts[[loc]] <- as.data.frame(table(location=rep(loc,n),Host=m$Host))

}
tests <- do.call(rbind,tests)
tests$p_BH_two_sided <- p.adjust(tests$p_two_sided,method="BH")
tests$BH_family <- "three_locations_times_two_viral_datasets"
tests$significant_BH_0_05 <- tests$p_BH_two_sided<0.05
tests$symbol <- ifelse(tests$p_BH_two_sided<0.001,"***",
  ifelse(tests$p_BH_two_sided<0.01,"**",ifelse(tests$p_BH_two_sided<0.05,"*","ns")))
write.csv(tests,file.path(out,"within_location_tests.csv"),row.names=FALSE,na="")
write.csv(do.call(rbind,summaries),file.path(out,"within_location_summary.csv"),row.names=FALSE)
write.csv(do.call(rbind,host_counts),file.path(out,"host_by_location.csv"),row.names=FALSE)
````
<!-- CODE within-location-vc-sharing-pair-comparisons END -->

<a id="vc-sharing-permutation-tests-for-all-samples-and-overlapping-locations"></a>

### Three-category VC-sharing tests for 69 samples

<!-- CODE vc-sharing-permutation-tests-for-all-samples-and-overlapping-locations BEGIN -->
````r
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop("Usage: Rscript script.R metadata_69.csv sample_pairwise_Jaccard.csv out_dir")
metadata_file <- args[1]
all_pair_file <- args[2]
out_dir <- args[3]

meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
all_pairs <- read.csv(all_pair_file, stringsAsFactors = FALSE)
datasets <- c("corrected_phage", "corrected_other")
stopifnot(nrow(meta)==69L,!anyDuplicated(meta$sample),!anyNA(meta),
  setequal(meta$location,c("JL_DH","YN_HH","YN_JH","YN_NJ","YN_PE","YN_YX")),
  nrow(all_pairs)==4692L,setequal(unique(all_pairs$dataset),datasets))


subgenus <- c(Aand = "Micrapis", Aflo = "Micrapis",
              Acer = "Apis", Amel = "Apis",
              Ador = "Megapis", Alab = "Megapis")

get_distances <- function(data, samples, index) {
  value <- sapply(datasets, function(dataset) {
    d <- data[data$dataset == dataset, ]
    mat <- matrix(NA_real_, length(samples), length(samples))
    i <- match(d$sample1, samples)
    j <- match(d$sample2, samples)
    mat[cbind(i, j)] <- d$Jaccard_distance
    mat[cbind(j, i)] <- d$Jaccard_distance
    mat[index]
  })
  stopifnot(all(is.finite(value)))
  value
}


two_sided_p <- function(observed, null, exact = FALSE) {
  upper <- sum(null >= observed - 1e-12)
  lower <- sum(null <= observed + 1e-12)
  if (exact) {
    min(1, 2 * min(upper, lower) / length(null))
  } else {
    min(1, 2 * (min(upper, lower) + 1) / (length(null) + 1))
  }
}

index69 <- t(combn(seq_len(nrow(meta)), 2))
distance69 <- get_distances(all_pairs, meta$sample, index69)

pair_category <- function(host, index) {
  left <- host[index[, 1]]
  right <- host[index[, 2]]
  ifelse(left == right, 1L,
         ifelse(subgenus[left] == subgenus[right], 2L, 3L))
}
category_names <- c("Same host", "Different hosts within subgenus",
                    "Cross-subgenus")
category69 <- pair_category(meta$Host, index69)
contrasts <- rbind(c(1L, 2L), c(1L, 3L), c(2L, 3L))

category_means <- function(category) {
  t(vapply(1:3, function(g) {
    colMeans(distance69[category == g, , drop = FALSE])
  }, numeric(length(datasets))))
}
means69 <- category_means(category69)


set.seed(20261025L)
n_perm <- 9999L
location_indices <- split(seq_len(nrow(meta)), meta$location)
permutation69 <- t(replicate(n_perm, {
  index <- seq_len(nrow(meta))
  for (g in location_indices) {
    index[g] <- g[sample.int(length(g))]
  }
  index
}))

null69 <- array(NA_real_, c(nrow(permutation69), 3L, length(datasets)))
for (b in seq_len(nrow(permutation69))) {
  host_b <- meta$Host[permutation69[b, ]]
  mean_b <- category_means(pair_category(host_b, index69))
  null69[b, , ] <- mean_b[contrasts[, 2], ] - mean_b[contrasts[, 1], ]
}

results69 <- list()
for (k in seq_along(datasets)) {
  for (j in seq_len(nrow(contrasts))) {
    g1 <- contrasts[j, 1]
    g2 <- contrasts[j, 2]
    observed <- means69[g2, k] - means69[g1, k]
    results69[[length(results69) + 1L]] <- data.frame(
      dataset = datasets[k], n_samples = nrow(meta),
      category1 = category_names[g1], category2 = category_names[g2],
      n_pairs_category1 = sum(category69 == g1),
      n_pairs_category2 = sum(category69 == g2),
      mean_distance_category1 = means69[g1, k],
      mean_distance_category2 = means69[g2, k],
      mean_sharing_category1 = 1 - means69[g1, k],
      mean_sharing_category2 = 1 - means69[g2, k],
      sharing_difference = observed,
      p_two_sided = two_sided_p(observed, null69[, j, k]),
      permutations = nrow(permutation69), permutation_block = "location")
  }
}
results69 <- do.call(rbind, results69)
results69$p_BH_two_sided <- p.adjust(results69$p_two_sided, method = "BH")

dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
write.csv(results69,file.path(out_dir,"analysis1_69_samples_three_categories.csv"),row.names=FALSE)
````
<!-- CODE vc-sharing-permutation-tests-for-all-samples-and-overlapping-locations END -->

<a id="module-03"></a>

## 03 Diversity analyses and phage-bacterium associations

<a id="audit-submission-abundance-reference-and-records"></a>


### Bacterial-phage Mantel tests and shared permutations

<!-- CODE bacteria-phage-community-association-sensitivity BEGIN -->
````r
suppressPackageStartupMessages({
  library(readxl)
  library(permute)
})
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[1],mustWork=TRUE) else normalizePath(".",mustWork=TRUE)
out <- if(length(args)>1L) args[2] else file.path(root,"mantel_results")
dir.create(out,showWarnings=FALSE,recursive=TRUE)

samples <- trimws(readLines(file.path(root,"samples.txt"),warn=FALSE))
samples <- samples[nzchar(samples)]
meta <- read.csv(file.path(root,"metadata_69.csv"),stringsAsFactors=FALSE)
stopifnot(identical(names(meta),c("sample","Host","location")),!anyDuplicated(meta$sample),
  !anyNA(meta),setequal(meta$location,c("JL_DH","YN_HH","YN_JH","YN_NJ","YN_PE","YN_YX")))
rownames(meta) <- meta$sample

bacteria_raw <- as.data.frame(read_excel(file.path(root,"correlation_new","Apis_bacteria_cvg_pct.xlsx")))
rownames(bacteria_raw) <- bacteria_raw[[1]]
bacteria <- t(as.matrix(bacteria_raw[,-1,drop=FALSE]))
storage.mode(bacteria) <- "double"

if(length(args)<3L) stop("Supply the phage matrix explicitly; see INPUTS.md for historical/new-matrix lineage.")
phage_raw <- read.csv(args[3],
  row.names=1,check.names=FALSE)
phage <- t(as.matrix(phage_raw))
storage.mode(phage) <- "double"

common <- samples[samples %in% rownames(meta) & samples %in% rownames(bacteria) &
  samples %in% rownames(phage)]
stopifnot(length(samples)==69,length(common)==69,!anyDuplicated(common))
meta <- meta[common,,drop=FALSE]
bacteria <- bacteria[common,,drop=FALSE]
phage <- phage[common,,drop=FALSE]

aitchison <- function(x,pseudocount=1e-6) {
  lx <- log(x+pseudocount)
  dist(sweep(lx,1,rowMeans(lx),"-"))
}
Db <- aitchison(bacteria)
Dp <- aitchison(phage)
pb <- cmdscale(Db,k=2,eig=TRUE)
pp <- cmdscale(Dp,k=2,eig=TRUE)
positive_b <- pb$eig[pb$eig>0]
positive_p <- pp$eig[pp$eig>0]

set.seed(20261024)
nperm <- 9999L
p_unrestricted <- shuffleSet(length(common),nset=nperm,control=how(nperm=nperm))
host_control <- how(nperm=nperm)
setBlocks(host_control) <- factor(meta$Host)
p_within_host <- shuffleSet(length(common),nset=nperm,control=host_control)
permutations <- list(unrestricted=p_unrestricted,within_Host=p_within_host)

average_rank <- function(v) {
  o <- order(v,na.last=NA)
  z <- v[o]
  groups <- cumsum(c(TRUE,diff(z)!=0))
  ranks <- ave(seq_along(z),groups,FUN=mean)
  out <- numeric(length(v));out[o] <- ranks;out
}
lower <- lower.tri(as.matrix(Db))
rank_b <- average_rank(as.matrix(Db)[lower])
rank_p <- average_rank(as.matrix(Dp)[lower])
rank_b <- rank_b-mean(rank_b)
rank_p <- rank_p-mean(rank_p)
rank_p_matrix <- matrix(0,length(common),length(common))
rank_p_matrix[lower] <- rank_p
rank_p_matrix <- rank_p_matrix+t(rank_p_matrix)
mantel_observed <- sum(rank_b*rank_p)/sqrt(sum(rank_b^2)*sum(rank_p^2))

results <- lapply(names(permutations),function(scheme) {
  p <- permutations[[scheme]]
  mantel_perm <- vapply(seq_len(nrow(p)),function(i) {
    y <- rank_p_matrix[p[i,],p[i,]][lower]
    sum(rank_b*y)/sqrt(sum(rank_b^2)*sum(y^2))
  },numeric(1))
  data.frame(analysis="Mantel_Aitchison_distance",scheme=scheme,
    statistic=mantel_observed,p=(1+sum(mantel_perm>=mantel_observed))/(1+nrow(p)),
    permutations=nrow(p))
})
results <- do.call(rbind,results)

variance <- rbind(
  data.frame(community="bacteria",axis=c("PCoA1","PCoA2"),
    variance_percent=100*pb$eig[1:2]/sum(positive_b)),
  data.frame(community="phage",axis=c("PCoA1","PCoA2"),
    variance_percent=100*pp$eig[1:2]/sum(positive_p))
)
scores <- data.frame(sample=common,Host=meta$Host,location=meta$location,
  bacteria_PCoA1=pb$points[,1],bacteria_PCoA2=pb$points[,2],
  phage_PCoA1=pp$points[,1],phage_PCoA2=pp$points[,2],check.names=FALSE)
write.csv(results,file.path(out,"association_sensitivity_results.csv"),row.names=FALSE)
write.csv(variance,file.path(out,"pcoa_variance_explained.csv"),row.names=FALSE)
write.csv(scores,file.path(out,"sample_PCoA12_scores.csv"),row.names=FALSE)
saveRDS(list(unrestricted=p_unrestricted,within_Host=p_within_host),
  file.path(out,"permutation_indices.rds"))
````
<!-- CODE bacteria-phage-community-association-sensitivity END -->

<a id="full-axis-bacterial-phage-procrustes"></a>

### Full-axis Procrustes test of bacterial and phage communities

<!-- CODE full-axis-bacterial-phage-procrustes BEGIN -->
````r
suppressPackageStartupMessages({library(readxl); library(vegan)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=6L) stop("Usage: Rscript script.R samples.txt bacteria.xlsx phage.csv metadata.csv permutations.rds out_dir")
inputs <- args[1:5]
out <- args[6]
stopifnot(all(file.exists(inputs)))
dir.create(out, recursive=TRUE, showWarnings=FALSE)
samples <- trimws(readLines(inputs[1], warn=FALSE))
samples <- samples[nzchar(samples)]
stopifnot(length(samples)==69L, !anyDuplicated(samples))

raw_b <- as.data.frame(read_excel(inputs[2]))
stopifnot(!anyDuplicated(raw_b[[1]]))
rownames(raw_b) <- raw_b[[1]]
bacteria <- t(as.matrix(raw_b[, -1, drop=FALSE]))
phage <- t(as.matrix(read.csv(inputs[3], row.names=1, check.names=FALSE)))
storage.mode(bacteria) <- storage.mode(phage) <- "double"
meta <- read.csv(inputs[4], stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(meta$sample), all(samples %in% meta$sample),
  all(samples %in% rownames(bacteria)), all(samples %in% rownames(phage)))
meta <- meta[match(samples, meta$sample), , drop=FALSE]
bacteria <- bacteria[samples, , drop=FALSE]
phage <- phage[samples, , drop=FALSE]
stopifnot(all(is.finite(bacteria)), all(is.finite(phage)),
  all(bacteria>=0), all(phage>=0), !anyNA(meta$Host),
  identical(names(meta),c("sample","Host","location")),
  setequal(meta$location,c("JL_DH","YN_HH","YN_JH","YN_NJ","YN_PE","YN_YX")))
permutations <- readRDS(inputs[5])
stopifnot(all(c("unrestricted", "within_Host") %in% names(permutations)))
permutations <- permutations[c("unrestricted", "within_Host")]
for(scheme in names(permutations)) {
  p <- permutations[[scheme]]
  stopifnot(is.matrix(p),nrow(p)==9999L,ncol(p)==69L,
    all(apply(p,1,function(z) identical(sort(as.integer(z)),seq_len(69L)))))
  if(scheme=="within_Host") stopifnot(all(apply(p,1,function(z) all(meta$Host[z]==meta$Host))))
}

aitchison <- function(x) {
  lx <- log(x + 1e-6)
  dist(sweep(lx, 1, rowMeans(lx), "-"))
}
Db <- aitchison(bacteria)
Dp <- aitchison(phage)
full_pcoa <- function(D, community) {
  d <- as.matrix(D)
  n <- nrow(d)
  J <- diag(n) - matrix(1/n, n, n)
  ev <- eigen(-0.5 * J %*% (d^2) %*% J, symmetric=TRUE)
  tolerance <- max(ev$values) * n * .Machine$double.eps * 10
  keep <- ev$values > tolerance
  stopifnot(!any(ev$values < -tolerance), any(keep))
  points <- sweep(ev$vectors[, keep, drop=FALSE], 2,
    sqrt(ev$values[keep]), "*")
  rownames(points) <- attr(D, "Labels")
  colnames(points) <- paste0("PCoA", seq_len(ncol(points)))
  error <- max(abs(as.matrix(dist(points)) - d))
  stopifnot(error / max(d) < 1e-8)
  retained <- sum(ev$values[keep]) / sum(ev$values[ev$values > 0])
  list(points=points, audit=data.frame(community=community, samples=n,
    effective_axes=sum(keep), retained_variance_percent=100*retained,
    eigenvalue_tolerance=tolerance, min_eigenvalue=min(ev$values),
    max_distance_reconstruction_error=error,
    relative_distance_reconstruction_error=error/max(d)),
    eigenvalues=NULL)
}
pb <- full_pcoa(Db, "bacteria")
pp <- full_pcoa(Dp, "phage")

fits <- lapply(names(permutations), function(scheme) {
  fit <- protest(pb$points, pp$points, permutations=permutations[[scheme]])
  manual_p <- (1+sum(fit$t >= fit$t0-sqrt(.Machine$double.eps))) /
    (1+length(fit$t))
  stopifnot(abs(fit$signif-manual_p)<1e-12)
  list(fit=fit, result=data.frame(analysis="Procrustes_all_effective_PCoA_axes",
    scheme=scheme, statistic=as.numeric(fit$t0), p=fit$signif,
    permutations=length(fit$t), m2=fit$ss,
    bacteria_axes=ncol(pb$points), phage_axes=ncol(pp$points),
    bacteria_variance_percent=pb$audit$retained_variance_percent,
    phage_variance_percent=pp$audit$retained_variance_percent))
})
names(fits) <- names(permutations)
results <- do.call(rbind, lapply(fits, `[[`, "result"))
write.csv(results,file.path(out,"full_axis_procrustes_results.csv"),row.names=FALSE)
saveRDS(list(distance_bacteria=Db,distance_phage=Dp,
  coordinates_bacteria=pb$points,coordinates_phage=pp$points,
  tests=lapply(fits,`[[`,"fit"),samples=samples,
  permutation_source=inputs[5],pseudocount=1e-6),file.path(out,"full_axis_analysis.rds"))
````
<!-- CODE full-axis-bacterial-phage-procrustes END -->


<a id="within-species-full-axis-bacterial-phage-procrustes"></a>

### Full-axis bacterial-phage Procrustes separately within six host species

<!-- CODE within-species-full-axis-bacterial-phage-procrustes BEGIN -->
````r
suppressPackageStartupMessages(library(vegan))
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=3L) stop("Usage: Rscript script.R full_axis_analysis.rds metadata_69.csv out_dir")
out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source <- readRDS(args[1])
meta <- read.csv(args[2],stringsAsFactors=FALSE)
meta <- meta[match(source$samples, meta$sample), c("sample", "Host", "location")]
hosts <- c("Aand", "Aflo", "Acer", "Amel", "Ador", "Alab")
species <- c("A. andreniformis", "A. florea", "A. cerana", "A. mellifera",
             "A. dorsata", "A. laboriosa")
stopifnot(nrow(meta) == 69, !anyNA(meta), !anyDuplicated(meta$sample),
          all(meta$Host %in% hosts), all(meta$location %in%
            c("JL_DH", "YN_HH", "YN_JH", "YN_NJ", "YN_PE", "YN_YX")))

# Subset the original distances, then reconstruct the complete ordination separately.
# Retain the same numerical eigenvalue tolerance as the original full-axis analysis.
full_pcoa <- function(D) {
  d <- as.matrix(D)
  n <- nrow(d)
  J <- diag(n) - matrix(1 / n, n, n)
  ev <- eigen(-.5 * J %*% (d^2) %*% J, symmetric = TRUE)
  tolerance <- max(ev$values) * n * .Machine$double.eps * 10
  keep <- ev$values > tolerance
  stopifnot(any(keep), !any(ev$values < -tolerance))
  coordinates <- sweep(ev$vectors[, keep, drop = FALSE], 2,
                        sqrt(ev$values[keep]), "*")
  rownames(coordinates) <- rownames(d)
  colnames(coordinates) <- paste0("PCoA", seq_len(ncol(coordinates)))
  error <- max(abs(as.matrix(dist(coordinates)) - d)) / max(d)
  stopifnot(error < 1e-8)
  list(points = coordinates,
       variance_percent = 100 * sum(ev$values[keep]) / sum(ev$values[ev$values > 0]),
       eigenvalue_tolerance = tolerance, relative_distance_error = error)
}

Db <- as.matrix(source$distance_bacteria)
Dp <- as.matrix(source$distance_phage)
fits <- results <- residuals <- list()
for (i in seq_along(hosts)) {
  host <- hosts[i]
  ids <- meta$sample[meta$Host == host]
  n <- length(ids)
  bacteria <- full_pcoa(as.dist(Db[ids, ids]))
  phages <- full_pcoa(as.dist(Dp[ids, ids]))

  # Both matrices contain only this host's samples, so every swap stays within host.
  # Use a separate reproducible seed for each host and 9,999 matched-sample swaps.
  seed <- 20261109L + i
  set.seed(seed)
  permutations <- t(replicate(9999L, sample.int(n)))
  fit <- protest(bacteria$points, phages$points, permutations = permutations)
  manual_p <- (1 + sum(fit$t >= fit$t0 - sqrt(.Machine$double.eps))) /
               (1 + length(fit$t))
  stopifnot(abs(manual_p - fit$signif) < 1e-12)
  results[[host]] <- data.frame(Host = host, species = species[i], samples = n,
    bacteria_axes = ncol(bacteria$points), phage_axes = ncol(phages$points),
    bacteria_variance_percent = bacteria$variance_percent,
    phage_variance_percent = phages$variance_percent,
    correlation = as.numeric(fit$t0), p_raw = fit$signif,
    permutations = length(fit$t), seed = seed,
    bacteria_relative_distance_error = bacteria$relative_distance_error,
    phage_relative_distance_error = phages$relative_distance_error)
  sq <- rowSums((fit$X-fit$Yrot)^2)
  stopifnot(length(sq)==length(ids),abs(sum(sq)-fit$ss)<1e-10)
  residuals[[host]] <- data.frame(sample=ids,Host=host,residual=sqrt(sq))
  fits[[host]] <- list(samples = ids, bacteria = bacteria, phages = phages,
                       fit = fit, permutations = permutations, seed = seed)
}
results <- do.call(rbind, results)
rownames(results) <- NULL
results$p_BH <- p.adjust(results$p_raw, method = "BH")
results$significant_BH <- results$p_BH < .05
write.csv(results, file.path(out, "Procrustes_by_host_results.csv"), row.names = FALSE)
write.csv(do.call(rbind,residuals),file.path(out,"Procrustes_by_host_residuals.csv"),row.names=FALSE)
write.csv(meta, file.path(out, "samples_used.csv"), row.names = FALSE)
saveRDS(fits, file.path(out, "Procrustes_by_host_models.rds"))
````
<!-- CODE within-species-full-axis-bacterial-phage-procrustes END -->

<a id="core-lefse-analysis"></a>

### Core LEfSe analysis

<!-- CODE core-lefse-analysis BEGIN -->
````bash
lefse_format_input.py abundance.lefse.txt host.in -c 2 -s 3 -u 1 -o 1000000
lefse_run.py host.in host.LDA2.res -l 2
````
<!-- CODE core-lefse-analysis END -->

<a id="diversity-mantel-tests-and-significance-letters"></a>

### Alpha diversity, significance letters and correlations

<!-- CODE diversity-mantel-tests-and-significance-letters BEGIN -->
````r
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("Usage: Rscript script.R input_bundle.rds output.rds")
bundle <- readRDS(args[1])
required <- c("meta","bacteria","phage","virus","eukvirus")
stopifnot(is.list(bundle), all(required %in% names(bundle)))
list2env(bundle[required], envir=environment())
stopifnot(identical(names(meta), c("sample","Host","location")),
  nrow(meta)==69L, !anyNA(meta), !anyDuplicated(meta$sample),
  setequal(meta$location,c("JL_DH","YN_HH","YN_JH","YN_NJ","YN_PE","YN_YX")))
samples <- meta$sample
for(M in list(bacteria,phage,virus,eukvirus)) stopifnot(
  !anyDuplicated(rownames(M)),all(samples %in% rownames(M)),all(is.finite(M)),all(M>=0))
bacteria <- bacteria[samples,,drop=FALSE]
phage <- phage[samples,,drop=FALSE]
virus <- virus[samples,,drop=FALSE]
eukvirus <- eukvirus[samples,,drop=FALSE]

calc_alpha <- function(M) data.frame(
  sample=rownames(M),shannon=vegan::diversity(M,index="shannon"),
  richness=rowSums(M>0))
alpha <- lapply(list(bacteria=bacteria,phage=phage,virus=virus,eukvirus=eukvirus),calc_alpha)
alpha_stats <- lapply(alpha,function(x) {
  x <- dplyr::left_join(x,meta,by="sample")
  list(kruskal=kruskal.test(shannon~Host,data=x),
       dunn=rstatix::dunn_test(x,shannon~Host,p.adjust.method="BH"))
})
dunn_p <- alpha_stats$phage$dunn
adjusted_p <- dunn_p$p.adj
names(adjusted_p) <- paste(dunn_p$group1,dunn_p$group2,sep="-")
significance_letters <- multcompView::multcompLetters(adjusted_p,threshold=0.05)$Letters

alpha_correlations <- list(
  pearson=cor.test(alpha$bacteria$shannon,alpha$phage$shannon,method="pearson"),
  spearman=cor.test(alpha$bacteria$shannon,alpha$phage$shannon,method="spearman"))
aitchison <- function(M) dist(compositions::clr(M+1e-6),method="euclidean")
D_phage <- as.matrix(aitchison(phage))
mean_within_host <- vapply(seq_along(samples),function(i) {
  same <- which(meta$Host==meta$Host[i])
  if(length(same)<2L) return(NA_real_)
  sum(D_phage[i,same])/(length(same)-1L)
},numeric(1))
alpha_beta_correlations <- list(
  spearman=cor.test(alpha$bacteria$shannon,mean_within_host,method="spearman"),
  pearson=cor.test(alpha$bacteria$shannon,mean_within_host,method="pearson"))
saveRDS(list(meta=meta,alpha=alpha,alpha_stats=alpha_stats,
  significance_letters=significance_letters,alpha_correlations=alpha_correlations,
  within_host_phage_distance=data.frame(sample=samples,mean_distance=mean_within_host),
  alpha_beta_correlations=alpha_beta_correlations),args[2])
````
<!-- CODE diversity-mantel-tests-and-significance-letters END -->

<a id="host-geography-marginal-permanova-and-dispersion"></a>

### Shared Aitchison distance preparation

<!-- CODE host-geography-marginal-permanova-and-dispersion BEGIN -->
````r
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else "."
out <- if(length(args)>1L) args[2] else file.path(root,"community_distances")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
samples <- scan(file.path(root,"samples.txt"),what=character(),quiet=TRUE)
meta <- read.csv(file.path(root,"metadata_69.csv"),stringsAsFactors=FALSE)
stopifnot(length(samples)==69L,!anyDuplicated(samples),!anyDuplicated(meta$sample),
          setequal(meta$sample,samples),!anyNA(meta))
meta <- meta[match(samples,meta$sample),,drop=FALSE]
for(nm in c("Host","location")) meta[[nm]] <- factor(meta[[nm]])
stopifnot(nlevels(meta$Host)==6L,setequal(meta$location,c("JL_DH","YN_HH","YN_JH","YN_NJ","YN_PE","YN_YX")))
rownames(meta) <- meta$sample
write_out <- function(x,name) write.csv(x,file.path(out,name),row.names=FALSE,na="")
write_out(meta,"metadata_69.csv")
write_out(as.data.frame(table(meta$Host,meta$location)),"host_by_location.csv")
if(length(args)!=5L) stop("Usage: Rscript script.R project_root out_dir all_virus.csv phage.csv other.csv")
paths <- setNames(args[3:5],c("all_viruses","phages","other_viruses"))
distances <- list(); manifests <- list()
for(dataset in names(paths)) {
  dat <- read.csv(paths[[dataset]],check.names=FALSE)
  stopifnot(!anyDuplicated(dat[[1]]),setequal(names(dat)[-1],samples))
  x <- t(as.matrix(dat[,samples,drop=FALSE])); storage.mode(x) <- "double"
  stopifnot(all(is.finite(x)),all(x>=0),max(abs(rowSums(x)-100))<1e-5)
  clr <- sweep(log(x+1e-6),1,rowMeans(log(x+1e-6)),"-")
  distances[[dataset]] <- dist(clr)
  saveRDS(distances[[dataset]],file.path(out,paste0(dataset,"_Aitchison.rds")))
  manifests[[dataset]] <- data.frame(dataset=dataset,path=paths[[dataset]],samples=nrow(x),
    features=ncol(x),pseudocount=1e-6,abundance_scale="existing_percent_0_to_100")
}
write_out(do.call(rbind,manifests),"input_manifest.csv")
````
<!-- CODE host-geography-marginal-permanova-and-dispersion END -->


<a id="subset-pcoa-and-permanova-from-existing-aitchison-distances"></a>

### PCoA and PERMANOVA for 69 and 61 samples

<!-- CODE subset-pcoa-and-permanova-from-existing-aitchison-distances BEGIN -->
````r
suppressPackageStartupMessages({library(vegan); library(permute)})
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=3L) stop("Usage: Rscript script.R project_root out_dir distance_dir")
root <- args[1]
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# 使用现有69样本图对应的距离，只从三列metadata读取地点。
meta <- read.csv(file.path(root, "metadata_69.csv"), stringsAsFactors = FALSE)
samples <- scan(file.path(root, "samples.txt"), what = "", quiet = TRUE)
meta <- meta[match(samples, meta$sample), c("sample", "Host", "location")]
locations <- c("JL_DH", "YN_HH", "YN_JH", "YN_NJ", "YN_PE", "YN_YX")
stopifnot(nrow(meta) == 69L, !anyNA(meta), !anyDuplicated(meta$sample),
          setequal(meta$location, locations))
files <- c(all_viruses = "all_viruses_Aitchison.rds", phages = "phages_Aitchison.rds",
           other_viruses = "other_viruses_Aitchison.rds")
full_distances <- lapply(files, function(f) {
  distance_dir <- if(length(args)>2L) args[3] else file.path(root,"host_geography_permanova_20261023")
  d <- readRDS(file.path(distance_dir,f))
  stopifnot(setequal(attr(d, "Labels"), meta$sample))
  d
})

# 子集只改变样本范围，保留原距离的完整特征空间。
scopes <- list(all_69=rep(TRUE,nrow(meta)),Yunnan_61=grepl("^YN_",meta$location))
stopifnot(sum(scopes$Yunnan_61)==61L)
results <- scores <- variance <- counts <- list()
for (si in seq_along(scopes)) {
  scope <- names(scopes)[si]
  m <- meta[scopes[[si]], ]
  m$Host <- factor(m$Host)
  m$location <- factor(m$location)
  counts[[scope]] <- data.frame(scope = scope,
    as.data.frame(table(Host = m$Host, location = m$location)))

  # 普通置换检验各模型；地点内置换只用于检验Host | Location。
  set.seed(20261023L + si * 100L + 1L)
  unrestricted <- shuffleSet(nrow(m), nset = 9999L,
    control = how(nperm = 9999L, within = Within(type = "free")))
  set.seed(20261023L + si * 100L + 2L)
  within_location <- shuffleSet(nrow(m), nset = 9999L,
    control = how(nperm = 9999L, within = Within(type = "free"), blocks = m$location))
  for (dataset in names(full_distances)) {
    d <- as.dist(as.matrix(full_distances[[dataset]])[m$sample, m$sample])
    pc <- cmdscale(d, eig = TRUE, k = 2)
    scores[[paste(scope, dataset)]] <- data.frame(scope = scope, dataset = dataset,
      sample = m$sample, Host = as.character(m$Host), location = as.character(m$location),
      PCoA1 = pc$points[, 1], PCoA2 = pc$points[, 2], row.names = NULL)
    variance[[paste(scope, dataset)]] <- data.frame(scope = scope, dataset = dataset,
      axis = 1:2, variance_percent = 100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0]))

    fits <- list(
      marginal_unrestricted = adonis2(d ~ Host + location, data = m,
        permutations = unrestricted, by = "margin"),
      Host_within_location = adonis2(d ~ Host + location, data = m,
        permutations = within_location, by = "margin"),
      joint_unrestricted = adonis2(d ~ Host + location, data = m,
        permutations = unrestricted))
    terms <- list(marginal_unrestricted = c("Host", "location"),
                  Host_within_location = "Host", joint_unrestricted = "Model")
    for (test in names(fits)) {
      tab <- as.data.frame(fits[[test]])[terms[[test]], , drop = FALSE]
      results[[paste(scope, dataset, test)]] <- data.frame(scope = scope,
        dataset = dataset, test = test, term = rownames(tab), samples = nrow(m),
        observed_locations = nlevels(m$location), permutations = 9999L,
        Df = tab$Df, R2 = tab$R2, pseudo_F = tab$F, p = tab[["Pr(>F)"]], row.names = NULL)
    }
  }
}
write.csv(do.call(rbind, results), file.path(out, "PERMANOVA_results.csv"), row.names = FALSE)
write.csv(do.call(rbind, scores), file.path(out, "PCoA_scores.csv"), row.names = FALSE)
write.csv(do.call(rbind, variance), file.path(out, "PCoA_variance.csv"), row.names = FALSE)
write.csv(do.call(rbind, counts), file.path(out, "host_by_location_counts.csv"), row.names = FALSE)
write.csv(meta, file.path(out, "metadata_69.csv"), row.names = FALSE)
````
<!-- CODE subset-pcoa-and-permanova-from-existing-aitchison-distances END -->

<a id="per-location-host-permanova-and-dispersion"></a>

### Host PERMANOVA within three locations

<!-- CODE per-location-host-permanova-and-dispersion BEGIN -->
````r
suppressPackageStartupMessages(library(vegan))
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=3L) stop("Usage: Rscript script.R project_root out_dir distance_dir")
root <- args[1]
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# 目的：在各多宿主地点内单独比较Host，沿用原图的距离与样本范围。
meta <- read.csv(file.path(root, "metadata_69.csv"), stringsAsFactors = FALSE)
samples <- scan(file.path(root, "samples.txt"), what = "", quiet = TRUE)
meta <- meta[match(samples, meta$sample), c("sample", "Host", "location")]
stopifnot(nrow(meta) == 69L, !anyNA(meta), !anyDuplicated(meta$sample),
  setequal(meta$location, c("JL_DH", "YN_HH", "YN_JH", "YN_NJ", "YN_PE", "YN_YX")))
files <- c(all_viruses = "all_viruses_Aitchison.rds", phages = "phages_Aitchison.rds",
           other_viruses = "other_viruses_Aitchison.rds")
full_distances <- lapply(files, function(f) {
  distance_dir <- if(length(args)>2L) args[3] else file.path(root,"host_geography_permanova_20261023")
  d <- readRDS(file.path(distance_dir,f))
  stopifnot(setequal(attr(d, "Labels"), meta$sample))
  d
})
locations <- c("JL_DH", "YN_JH", "YN_YX")
n_perm <- 9999L
host_tests <- counts_out <- models <- list()

for (li in seq_along(locations)) {
  loc <- locations[li]
  m <- meta[meta$location == loc, ]
  m$Host <- factor(m$Host)
  host <- as.integer(m$Host)
  n <- nrow(m)
  counts <- tabulate(host)
  counts_out[[loc]] <- data.frame(location = loc, Host = levels(m$Host), n_samples = counts)
  n_allocations <- round(exp(lgamma(n + 1) - sum(lgamma(counts + 1))))
  exact <- length(counts) == 2L && n_allocations <= n_perm
  set.seed(2034694L + li)

  if (exact) {
    # 两宿主的小地点：枚举第一宿主占据的全部位置，每种标签分配保留一次。
    chosen <- t(combn(seq_len(n), counts[1]))
    allocations <- matrix(2L, nrow(chosen), n)
    allocations[cbind(rep(seq_len(nrow(chosen)), each = ncol(chosen)),
                       as.vector(t(chosen)))] <- 1L
    permutations <- t(vapply(seq_len(nrow(allocations)), function(b) {
      index <- integer(n)
      for (g in seq_along(counts)) index[host == g] <- which(allocations[b, ] == g)
      index
    }, integer(n)))
    # adonis2会将观测结果加回参考集，因此只提交其余标签分配。
    observed_row <- apply(permutations, 1, function(x) all(x == seq_len(n)))
    permutations <- permutations[!observed_row, , drop = FALSE]
    method <- "Exact host-label PERMANOVA"
  } else {
    permutations <- t(replicate(n_perm, sample.int(n)))
    method <- "Random-permutation PERMANOVA"
  }

  for (dataset in names(full_distances)) {
    # 地点内Location为常数；模型仅包含Host，不重新筛选距离的特征空间。
    d <- as.dist(as.matrix(full_distances[[dataset]])[m$sample, m$sample])
    fit <- adonis2(d ~ Host, data = m, permutations = permutations, by = "margin")
    tab <- as.data.frame(fit)["Host", , drop = FALSE]
    host_tests[[paste(loc, dataset)]] <- data.frame(location = loc, dataset = dataset,
      n_samples = n, n_hosts = nlevels(m$Host), Df = tab$Df,
      residual_Df = fit["Residual", "Df"], R2 = tab$R2, pseudo_F = tab$F,
      p = tab[["Pr(>F)"]], method = method, permutations_supplied = nrow(permutations),
      reference_size_including_observed = nrow(permutations) + 1L)

    pc <- cmdscale(d,eig=TRUE,k=2)
    models[[paste(loc,dataset)]] <- list(metadata=m,distance=d,host_test=fit,pcoa=pc)

  }
}

host_results <- do.call(rbind, host_tests)
host_results$p_BH_9_tests <- p.adjust(host_results[["p"]], method = "BH")
write.csv(host_results,file.path(out,"PERMANOVA_by_location.csv"),row.names=FALSE)
write.csv(do.call(rbind,counts_out),file.path(out,"sample_counts_by_host_location.csv"),row.names=FALSE)
saveRDS(models,file.path(out,"model_objects.rds"))
````
<!-- CODE per-location-host-permanova-and-dispersion END -->

<a id="path-models-and-effect-sizes"></a>

### Phage path models and effect sizes

<!-- CODE path-models-and-effect-sizes BEGIN -->
````r
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("Usage: Rscript script.R input_bundle.rds output.rds")
bundle <- readRDS(args[1])
required <- c("meta","bacteria","phage","alpha")
stopifnot(is.list(bundle), all(required %in% names(bundle)))
list2env(bundle[required], envir=environment())
stopifnot(identical(names(meta), c("sample","Host","location")),
  nrow(meta)==69L, !anyNA(meta), !anyDuplicated(meta$sample),
  setequal(meta$location,c("JL_DH","YN_HH","YN_JH","YN_NJ","YN_PE","YN_YX")))
aitchison <- function(M) dist(compositions::clr(M+1e-6),method="euclidean")
samples <- meta$sample
for(M in list(bacteria,phage)) stopifnot(!anyDuplicated(rownames(M)),
  all(samples %in% rownames(M)),all(is.finite(M)),all(M>=0))
bacteria <- bacteria[samples,,drop=FALSE]
phage <- phage[samples,,drop=FALSE]
stopifnot(identical(alpha$bacteria$sample,samples),identical(alpha$phage$sample,samples))
fit_paths <- function(bacterial_score,viral_score,meta) {
  dat <- transform(meta,
    z_bacteria=as.numeric(scale(bacterial_score)),
    z_virus=as.numeric(scale(viral_score)))
  dat$Host <- factor(dat$Host)
  dat$location <- factor(dat$location)

  # Keep both categorical predictors; use separate dummy variables, not ranks.
  X <- model.matrix(~Host+location,data=dat)[,-1,drop=FALSE]
  colnames(X) <- make.names(colnames(X))
  predictors <- colnames(X)
  observed <- data.frame(X,Bacteria=dat$z_bacteria,Virus=dat$z_virus)
  if(anyNA(observed) || qr(cbind(1,X))$rank!=ncol(X)+1L)
    stop("The complete-case path design must be finite and full rank")

  paths <- rbind(
    data.frame(from=predictors,to="Bacteria",parameter=paste0("b",seq_along(predictors))),
    data.frame(from=predictors,to="Virus",parameter=paste0("v",seq_along(predictors))),
    data.frame(from="Bacteria",to="Virus",parameter="bacteria_virus"))
  equations <- c(
    paste(paths$from,"->",paths$to,",",paths$parameter,", NA"),
    "Bacteria <-> Bacteria, residual_bacteria, NA",
    "Virus <-> Virus, residual_virus, NA")
  ram <- sem::specifyModel(text=paste(equations,collapse="\n"),quiet=TRUE)
  # Default objectiveML: fit the complete recursive observed-variable model.
  # Exogenous variances/covariances are fixed at their observed values.
  sem_model <- sem::sem(ram,S=cov(observed),N=nrow(observed),
                        fixed.x=predictors,robust=FALSE)
  if(!isTRUE(sem_model$convergence)) stop("SEM did not converge")
  estimate <- stats::coef(sem_model)
  standard_error <- sqrt(diag(stats::vcov(sem_model)))
  paths$estimate <- unname(estimate[paths$parameter])
  paths$SE <- unname(standard_error[paths$parameter])
  paths$z <- paths$estimate/paths$SE
  paths$p <- 2*stats::pnorm(-abs(paths$z))

  # Matching linear submodels retain the original Type II F tests and
  # effectsize partial-eta-squared calculation as supplementary statistics.
  # Their t/F p values are distinct from the SEM Wald-z p values above.
  m_bacteria <- stats::lm(z_bacteria~Host+location,data=dat)
  m_virus <- stats::lm(z_virus~z_bacteria+Host+location,data=dat)
  R2 <- c(bacteria=1-sem_model$P["Bacteria","Bacteria"]/var(dat$z_bacteria),
          virus=1-sem_model$P["Virus","Virus"]/var(dat$z_virus))
  adj_R2 <- 1-(1-R2)*(nrow(dat)-1)/
    c(stats::df.residual(m_bacteria),stats::df.residual(m_virus))
  main <- paths[paths$parameter=="bacteria_virus",]
  list(sem_model=sem_model,sem_paths=paths,
    standardized_path=c(Estimate=main$estimate,SE=main$SE,z=main$z,p=main$p),
    bacterial_model=m_bacteria,viral_model=m_virus,
    bacterial_typeII=car::Anova(m_bacteria,type=2),
    viral_typeII=car::Anova(m_virus,type=2),
    bacterial_partial_eta2=effectsize::eta_squared(m_bacteria,partial=TRUE),
    viral_partial_eta2=effectsize::eta_squared(m_virus,partial=TRUE),
    R2=R2,adj_R2=adj_R2,metadata=dat,
    global_model_df=stats::df.residual(sem_model))
}
sem_alpha_phage <- fit_paths(alpha$bacteria$shannon,alpha$phage$shannon,meta)
bacteria_PCoA1 <- cmdscale(aitchison(bacteria),k=2,eig=TRUE)$points[,1]
phage_PCoA1 <- cmdscale(aitchison(phage),k=2,eig=TRUE)$points[,1]
sem_beta_phage <- fit_paths(bacteria_PCoA1,phage_PCoA1,meta)
saveRDS(list(alpha=sem_alpha_phage,beta=sem_beta_phage),args[2])
````
<!-- CODE path-models-and-effect-sizes END -->


<a id="module-04"></a>

## 04 Bacterial host prediction for bacteriophages

DIAMOND/BLASTN host transfer reuses earlier alignment results; CHERRY host predictions are taken from the PhaBOX v2 web server output.

<a id="host-prediction-trna-matching"></a>

### Exact tRNA matching

<!-- CODE host-prediction-trna-matching BEGIN -->
````bash
set -euo pipefail
bacterial_contigs=$1
phage_contigs=$2
out=$3
mkdir -p "$out"

aragorn -t -w -fasta -o "$out/bacteria.tRNA.fa" "$bacterial_contigs"
aragorn -t -w -fasta -o "$out/phage.tRNA.fa" "$phage_contigs"
makeblastdb -in "$out/bacteria.tRNA.fa" -dbtype nucl -out "$out/bacteria.tRNA"
blastn \
    -db "$out/bacteria.tRNA" \
    -query "$out/phage.tRNA.fa" \
    -outfmt "6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -out "$out/phage_bacteria.tRNA.blastn.tsv" \
    -perc_identity 100 -qcov_hsp_perc 100 \
    -num_threads "${THREADS:-1}"
````
<!-- CODE host-prediction-trna-matching END -->

<a id="wish-models-and-background-likelihoods"></a>

### WIsH models and background likelihoods

<!-- CODE wish-models-and-background-likelihoods BEGIN -->
````bash
set -euo pipefail
bacterial_dir=$1
background_phage_dir=$2
model_dir=$3
null_dir=$4
mkdir -p "$model_dir" "$null_dir"

WIsH -c build -g "$bacterial_dir" -m "$model_dir" -t "${THREADS:-1}"
WIsH -c predict -g "$background_phage_dir" -m "$model_dir" \
    -r "$null_dir" -b 1 -t "${THREADS:-1}"
````
<!-- CODE wish-models-and-background-likelihoods END -->

<a id="wish-host-prediction-with-null-parameters"></a>

### WIsH host prediction

<!-- CODE wish-host-prediction-with-null-parameters BEGIN -->
````bash
set -euo pipefail
phage_dir=$1
model_dir=$2
null_parameters=$3
out=$4
test -s "$null_parameters" || { printf '%s\n' "Missing matching WIsH nullParameters.tsv" >&2; exit 1; }
mkdir -p "$out"

WIsH -c predict -g "$phage_dir" -m "$model_dir" -r "$out" \
    -n "$null_parameters" -b 1 -p -t "${THREADS:-1}"
````
<!-- CODE wish-host-prediction-with-null-parameters END -->

<a id="module-05"></a>

## 05 Prediction and functional annotation of viral proteins

<a id="prodigal-protein-prediction"></a>

### Prodigal protein prediction

<!-- CODE prodigal-protein-prediction BEGIN -->
````bash
prodigal -i final_virus.min1000.fa -a final_virus.min1000.faa \
    -d final_virus.min1000.cds -f gbk -o final_virus.min1000.gbk -g 11 -p meta
````
<!-- CODE prodigal-protein-prediction END -->

<a id="phrog-hh-suite-annotation"></a>

### PHROG HH-suite

<!-- CODE phrog-hh-suite-annotation BEGIN -->
````bash
seqkit replace -p ' .*$' -r '' final_virus.min1000.faa > proteins.clean.faa
ffindex_from_fasta -s proteins.ffdata proteins.ffindex proteins.clean.faa
hhblits_omp -cpu 1 -i proteins -d "$PHROG_DB" \
    -o proteins_phrogs -blasttab proteins_phrogs.tsv
tr -d '\000\r' < proteins_phrogs.tsv.ffdata > proteins_phrogs.tsv
````
<!-- CODE phrog-hh-suite-annotation END -->

<a id="phrog-best-hit-filtering"></a>

### PHROG best-hit filtering

<!-- CODE phrog-best-hit-filtering BEGIN -->
````bash
awk '
BEGIN {OFS="\t"}
/^#/ || NF < 12 {next}
{
  q  = $1     # query: 蛋白 ID
  ph = $2     # PHROG ID
  e  = $11 + 0  # E-value
  s  = $12 + 0  # bit score

  if (e > 1e-5) next       # E-value 阈值
  if (s < 40)  next        # Score 必须 >= 40

  if (!(q in best_e) || e < best_e[q] || (e == best_e[q] && s > best_s[q])) {
    best_e[q]  = e
    best_s[q]  = s
    best_ph[q] = ph
  }
}
END {
  OFS = "\t"
  print "protein_id", "PHROG_ID", "Evalue", "Score"
  for (q in best_ph) {
    print q, best_ph[q], best_e[q], best_s[q]
  }
}' proteins_phrogs.tsv \
> Apis_meta_allsample_with_prophage_min1000_phrogs.besthit.tsv
````
<!-- CODE phrog-best-hit-filtering END -->

<a id="gene-abundance-mapping"></a>

### Gene-abundance mapping

<!-- CODE gene-abundance-mapping BEGIN -->
````bash
while IFS= read -r sample; do
    2bwt-builder "${sample}.cds"
    soap -a "${sample}.clean.1.fq.gz" -b "${sample}.clean.2.fq.gz" \
        -D "${sample}.cds.index" -o "${sample}.genes.paired.aln" \
        -2 "${sample}.genes.single.aln" -M 4 -l 30 -r 1 -v 3 -m 200 -p 12
    soap.coverage -cvg -i "${sample}.genes.paired.aln" "${sample}.genes.single.aln" \
        -refsingle "${sample}.cds" -o "${sample}.gene_coverage.txt"
    gawk 'BEGIN{OFS="\t"}{
      gene=$1;sub(/:$/,"",gene);split($2,a,"/");len=a[1]+0
      depth=0;if(match($0,/Depth:([0-9.]+)/,m))depth=m[1]+0
      print gene,len,depth,len*depth
    }' "${sample}.gene_coverage.txt" > "${sample}.gene_abundance.tsv"
done < samples.txt
````
<!-- CODE gene-abundance-mapping END -->

<a id="phrog-functional-composition"></a>

### PHROG annotation rate

<!-- CODE phrog-functional-composition BEGIN -->
````r
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("Usage: Rscript script.R input_bundle.rds output.rds")
bundle <- readRDS(args[1])
required <- c("meta","phrog_gene_counts")
stopifnot(is.list(bundle), all(required %in% names(bundle)))
list2env(bundle[required], envir=environment())
stopifnot(identical(names(meta), c("sample","Host","location")),
  nrow(meta)==69L, !anyNA(meta), !anyDuplicated(meta$sample),
  setequal(meta$location,c("JL_DH","YN_HH","YN_JH","YN_NJ","YN_PE","YN_YX")))
phrog_annotation_rate <- phrog_gene_counts |>
  dplyr::mutate(annotation_rate=100*annotated_genes/total_genes)
saveRDS(phrog_annotation_rate,args[2])
````
<!-- CODE phrog-functional-composition END -->

<a id="module-06"></a>

## 06 Phylogenetic analysis

<a id="amfv-sequence-extraction-alignment-and-phylogeny"></a>

### AmFV sequence extraction, alignment, and phylogeny

<!-- CODE amfv-sequence-extraction-alignment-and-phylogeny BEGIN -->
````bash
for gene in Dda NrdA; do
    mafft --auto "${gene}.raw.faa" > "${gene}.aln.faa"
    iqtree2 -s "${gene}.aln.faa" -st AA -m LG+G4 \
        -bb 1000 -alrt 1000 -nt 8 -pre "$gene"
done
````
<!-- CODE amfv-sequence-extraction-alignment-and-phylogeny END -->


<a id="phrog-protein-phylogenies"></a>

### Protein phylogenies for the 18 reported PHROGs

<!-- CODE phrog-protein-phylogenies BEGIN -->
````bash
while IFS= read -r phrog; do
    mafft --auto "${phrog}.faa" > "${phrog}.aln.faa"
    iqtree2 -s "${phrog}.aln.faa" -st AA -m LG+G4 \
        -bb 1000 -alrt 1000 -nt 8 -pre "$phrog"
done < phrog_reported_ids.txt
````
<!-- CODE phrog-protein-phylogenies END -->

<a id="module-07"></a>

## 07 Phage-bacterium network construction and robustness analysis

<a id="core-network-filtering-topology-and-robustness-functions"></a>

### Network inference, topology and robustness

<!-- CODE core-network-filtering-topology-and-robustness-functions BEGIN -->
````r
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("Usage: Rscript script.R input_bundle.rds output.rds")
bundle <- readRDS(args[1])
required <- c("meta","bacteria","phage")
stopifnot(is.list(bundle), all(required %in% names(bundle)))
list2env(bundle[required], envir=environment())
stopifnot(identical(names(meta), c("sample","Host","location")),
  nrow(meta)==69L, !anyNA(meta), !anyDuplicated(meta$sample),
  setequal(meta$location,c("JL_DH","YN_HH","YN_JH","YN_NJ","YN_PE","YN_YX")))
samples <- meta$sample
for(M in list(bacteria,phage)) stopifnot(!anyDuplicated(rownames(M)),
  all(samples %in% rownames(M)),all(is.finite(M)),all(M>=0))
bacteria <- bacteria[samples,,drop=FALSE]
phage <- phage[samples,,drop=FALSE]
stopifnot(!any(colnames(bacteria) %in% colnames(phage)))
set.seed(1)
present_threshold <- 0.001
min_prevalence_prop <- 0.10
min_mean_rel_abund <- 0
max_attack_frac <- 0.8
random_trials <- 50L

filter_features <- function(M, max_nodes, host_tag) {
  rs <- rowSums(M, na.rm = TRUE)
  rs[rs == 0] <- 1

  R <- sweep(M, 1, rs, "/")

  prev <- colMeans(R > present_threshold, na.rm = TRUE)
  meanrel <- colMeans(R, na.rm = TRUE)

  keep <- which(prev >= min_prevalence_prop & meanrel >= min_mean_rel_abund)

  if (length(keep) == 0) {
    return(M[, FALSE, drop = FALSE])
  }

  Rk <- R[, keep, drop = FALSE]

  if (is.finite(max_nodes) && ncol(Rk) > max_nodes) {
    ord <- order(colMeans(Rk), decreasing = TRUE)
    Rk <- Rk[, ord[seq_len(max_nodes)], drop = FALSE]
  }

  Rk
}

prune_graph <- function(g, abs_min = 0) {
  if (igraph::ecount(g) == 0L) return(g)

  if (!"weight" %in% igraph::edge_attr_names(g)) {
    g <- igraph::set_edge_attr(g, "weight", value = rep(1, igraph::ecount(g)))
  }

  w <- abs(igraph::edge_attr(g, "weight"))

  drop_ids <- which(w < abs_min)
  if(length(drop_ids)) g <- igraph::delete_edges(g,drop_ids)

  g
}

drop_virus_virus_edges <- function(g, vir_ids) {
  if (igraph::ecount(g) == 0L || length(vir_ids) == 0L) return(g)

  Edf <- igraph::as_data_frame(g, what = "edges")
  v1_is_virus <- Edf[, 1] %in% vir_ids
  v2_is_virus <- Edf[, 2] %in% vir_ids

  keep <- !(v1_is_virus & v2_is_virus)

  if (all(keep)) return(g)

  igraph::graph_from_data_frame(
    Edf[keep, , drop = FALSE],
    directed = FALSE,
    vertices = data.frame(name = igraph::V(g)$name, stringsAsFactors = FALSE)
  )
}

net_stats <- function(g) {
  if (igraph::vcount(g) == 0 || igraph::ecount(g) == 0) {
    return(data.frame(
      n = igraph::vcount(g),
      m = igraph::ecount(g),
      density = 0,
      mean_deg = 0,
      sd_deg = 0,
      transitivity = NA_real_,
      assortativity = NA_real_,
      n_comp = if (igraph::vcount(g) == 0) 0 else igraph::components(g)$no,
      lcc_size = 0,
      lcc_frac = 0,
      avg_path_LCC = NA_real_,
      diameter_LCC = NA_real_,
      modularity = NA_real_,
      n_communities = NA_integer_,
      mean_abs_weight = NA_real_,
      frac_pos = NA_real_,
      frac_neg = NA_real_,
      kcore_max = NA_real_,
      kcore_mean = NA_real_,
      spectral_radius = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  deg <- igraph::degree(g)
  dens <- igraph::edge_density(g)
  tr <- igraph::transitivity(g, type = "global", isolates = "zero")

  asrt <- tryCatch(
    igraph::assortativity_degree(g, directed = FALSE),
    error = function(e) NA_real_
  )

  comps <- igraph::components(g)
  lcc_vids <- which(comps$membership == which.max(comps$csize))
  gL <- igraph::induced_subgraph(g, lcc_vids)

  apl <- tryCatch(
    igraph::average.path.length(gL, directed = FALSE),
    error = function(e) NA_real_
  )

  diam <- tryCatch(
    igraph::diameter(gL, directed = FALSE),
    error = function(e) NA_real_
  )

  gL_unw <- if ("weight" %in% igraph::edge_attr_names(gL)) {
    igraph::delete_edge_attr(gL, "weight")
  } else {
    gL
  }

  comm <- tryCatch(
    igraph::cluster_louvain(gL_unw),
    error = function(e) NULL
  )

  mod <- if (is.null(comm)) NA_real_ else igraph::modularity(comm)
  ncom <- if (is.null(comm)) NA_integer_ else length(unique(comm$membership))

  w <- igraph::edge_attr(g, "weight")
  mean_abs_w <- mean(abs(w))
  frac_pos <- mean(w > 0)
  frac_neg <- mean(w < 0)

  kc <- tryCatch(
    igraph::coreness(g),
    error = function(e) rep(NA_real_, igraph::vcount(g))
  )

  kcore_max <- suppressWarnings(max(kc, na.rm = TRUE))
  kcore_mean <- suppressWarnings(mean(kc, na.rm = TRUE))

  A <- as.matrix(igraph::as_adj(g, attr = "weight", sparse = FALSE))
  sr <- tryCatch(
    max(abs(eigen(A, symmetric = TRUE, only.values = TRUE)$values)),
    error = function(e) NA_real_
  )

  data.frame(
    n = igraph::vcount(g),
    m = igraph::ecount(g),
    density = dens,
    mean_deg = mean(deg),
    sd_deg = stats::sd(deg),
    transitivity = tr,
    assortativity = asrt,
    n_comp = comps$no,
    lcc_size = max(comps$csize),
    lcc_frac = max(comps$csize) / igraph::vcount(g),
    avg_path_LCC = apl,
    diameter_LCC = diam,
    modularity = mod,
    n_communities = ncom,
    mean_abs_weight = mean_abs_w,
    frac_pos = frac_pos,
    frac_neg = frac_neg,
    kcore_max = kcore_max,
    kcore_mean = kcore_mean,
    spectral_radius = sr,
    stringsAsFactors = FALSE
  )
}

lcc_fraction <- function(g) {
  if (igraph::vcount(g) == 0) return(0)
  cs <- igraph::components(g)$csize
  max(cs) / igraph::vcount(g)
}

attack_curve_lcc_single <- function(g0, strategy = c("random", "degree", "betweenness")) {
  strategy <- match.arg(strategy)

  g <- igraph::delete_vertices(g0, which(igraph::degree(g0) == 0))
  n0 <- igraph::vcount(g)

  if (n0 == 0) {
    return(data.frame(
      step = integer(0),
      frac_removed = numeric(0),
      lcc = numeric(0)
    ))
  }

  lccs <- numeric(n0 + 1)
  lccs[1] <- lcc_fraction(g)

  for (k in 1:n0) {
    if (strategy == "random") {
      v <- sample(igraph::V(g), 1)

    } else if (strategy == "degree") {
      d <- igraph::degree(g)
      vmax <- names(which(d == max(d)))[1]
      v <- igraph::V(g)[vmax]

    } else {
      b <- igraph::betweenness(g, directed = FALSE, weights = NA)
      vmax <- names(which(b == max(b)))[1]
      v <- igraph::V(g)[vmax]
    }

    g <- igraph::delete_vertices(g, v)
    lccs[k + 1] <- lcc_fraction(g)
  }

  data.frame(
    step = 0:n0,
    frac_removed = (0:n0) / n0,
    lcc = lccs
  )
}

attack_curve_lcc <- function(g, strategy = c("random", "degree", "betweenness"), trials = 50L) {
  strategy <- match.arg(strategy)

  if (strategy == "random") {
    acc <- NULL

    for (t in seq_len(trials)) {
      cc <- attack_curve_lcc_single(g, "random")
      acc <- if (is.null(acc)) cc else transform(acc, lcc = lcc + cc$lcc)
    }

    acc$lcc <- acc$lcc / trials
    acc$strategy <- "random"
    acc

  } else {
    cc <- attack_curve_lcc_single(g, strategy)
    cc$strategy <- strategy
    cc
  }
}

nc_bin_norm <- function(adj_matrix) {
  A <- as.matrix(adj_matrix)
  A[is.na(A)] <- 0
  A[abs(A) != 0] <- 1

  n <- nrow(A)
  if (n <= 1) return(0)

  ev <- eigen(A, only.values = TRUE)$values
  eta <- log(sum(exp(ev)) / n)

  K <- matrix(1, n, n)
  diag(K) <- 0
  evK <- eigen(K, only.values = TRUE)$values
  eta_max <- log(sum(exp(evK)) / n)

  as.numeric(eta / eta_max)
}

attack_curve_nnc <- function(
    g,
    strategy = c("random", "degree", "betweenness"),
    max_frac = max_attack_frac,
    random_trials = random_trials
) {
  strategy <- match.arg(strategy)

  g0 <- igraph::delete_vertices(g, which(igraph::degree(g) == 0))
  n0 <- igraph::vcount(g0)

  if (n0 == 0) {
    return(data.frame(
      step = integer(),
      frac_removed = numeric(),
      nnc = numeric()
    ))
  }

  max_rm <- floor(max_frac * n0)
  K <- if (strategy == "random") random_trials else 1L

  acc <- NULL

  for (rep in seq_len(K)) {
    g_cur <- g0

    steps <- 0:max_rm
    vals <- numeric(length(steps))
    vals[1] <- nc_bin_norm(igraph::as_adjacency_matrix(g_cur, sparse = FALSE))

    for (i in 1:max_rm) {
      if (igraph::vcount(g_cur) == 0) {
        vals[i + 1] <- 0
        next
      }

      if (strategy == "random") {
        v <- sample(igraph::V(g_cur), 1)

      } else if (strategy == "degree") {
        d <- igraph::degree(g_cur)
        v <- igraph::V(g_cur)[which.max(d)]

      } else {
        b <- igraph::betweenness(g_cur, directed = FALSE, weights = NA)
        v <- igraph::V(g_cur)[which.max(b)]
      }

      g_cur <- igraph::delete_vertices(g_cur, v)

      vals[i + 1] <- if (igraph::vcount(g_cur) == 0) {
        0
      } else {
        nc_bin_norm(igraph::as_adjacency_matrix(g_cur, sparse = FALSE))
      }
    }

    df <- data.frame(
      step = steps,
      frac_removed = steps / n0,
      nnc = vals
    )

    acc <- if (is.null(acc)) df else transform(acc, nnc = nnc + df$nnc)
  }

  if (K > 1) acc$nnc <- acc$nnc / K
  acc
}

infer_graph <- function(matrices) {
  fit <- SpiecEasi::spiec.easi(matrices,method="glasso",nlambda=300,
    lambda.min.ratio=0.3,
    pulsar.params=list(thresh=0.2,subsample.ratio=0.8,rep.num=50,seed=1))
  ids <- if(is.list(matrices))unlist(lapply(matrices,colnames))else colnames(matrices)
  W <- as.matrix(stats::cov2cor(as.matrix(SpiecEasi::getOptCov(fit)))*SpiecEasi::getRefit(fit))
  dimnames(W) <- list(ids,ids)
  g <- igraph::graph_from_adjacency_matrix(W,mode="undirected",diag=FALSE,weighted=TRUE)
  g <- prune_graph(g,abs_min=0.05)
  igraph::delete_vertices(g,which(igraph::degree(g)==0))
}
networks <- list()
for(host in c("Aand","Aflo","Acer","Amel","Ador","Alab")) {
  ids <- meta$sample[meta$Host==host]
  B <- filter_features(bacteria[ids,,drop=FALSE],120,paste0(host,"-bac"))
  P <- filter_features(phage[ids,,drop=FALSE],120,paste0(host,"-phage"))
  if(nrow(B)<3 || ncol(B)==0 || ncol(P)==0)next
  B <- B[,colSums(B)>0,drop=FALSE]
  P <- P[,colSums(P)>0,drop=FALSE]
  if(ncol(B)==0 || ncol(P)==0)next
  bac_graph <- infer_graph(B)
  bv_graph <- infer_graph(list(B,P))
  filtered_graph <- drop_virus_virus_edges(bv_graph,colnames(P))
  filtered_graph <- igraph::delete_vertices(filtered_graph,which(igraph::degree(filtered_graph)==0))
  if(igraph::ecount(bv_graph)>0) {
    community <- igraph::cluster_louvain(bv_graph,weights=abs(igraph::E(bv_graph)$weight))
    igraph::V(bv_graph)$community <- igraph::membership(community)
  }
  networks[[host]] <- list(bacteria=bac_graph,BV=bv_graph,BV_filtered=filtered_graph)
}
auc_rows <- list()
for(metric in c("LCC","NNC"))for(host in names(networks)) {
  for(net_name in c("bacteria","BV"))for(strategy in c("random","degree","betweenness")) {
    g <- networks[[host]][[net_name]]
    auc <- 0
    if(igraph::vcount(g)>=10) {
      curve <- if(metric=="LCC")attack_curve_lcc(g,strategy,trials=50L)else
        attack_curve_nnc(g,strategy,max_frac=0.8,random_trials=50L)
      curve <- curve[curve$frac_removed<=0.8,,drop=FALSE]
      y <- if(metric=="LCC")curve$lcc else curve$nnc
      if(nrow(curve)>=2)auc <- pracma::trapz(curve$frac_removed,y)
    }
    auc_rows[[length(auc_rows)+1]] <- data.frame(host,network=net_name,strategy,metric,AUC=auc)
  }
}
robustness_auc <- do.call(rbind,auc_rows)
topology <- list()
for(host in names(networks)) for(net_name in c("bacteria","BV")) {
  topology[[paste(host,net_name)]] <- cbind(host=host,network=net_name,
    net_stats(networks[[host]][[net_name]]))
}
saveRDS(list(networks=networks,topology=do.call(rbind,topology),robustness_auc=robustness_auc),args[2])
````
<!-- CODE core-network-filtering-topology-and-robustness-functions END -->
