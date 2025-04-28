#!/bin/bash

# 自动识别样本ID列表
samples=($(ls 00_raw_data/*_1.fastq.gz | xargs -n 1 basename | sed 's/_1.fastq.gz//' | sort | uniq))

# 线程数
threads=18

# 参考基因组和注释文件路径
star_index="/home/future/index/star_index_hg38"
rsem_ref="/home/future/index/rsem_reference_hg38/hg38_rsem"
gtf="/home/future/GRCh38_reference/Homo_sapiens.GRCh38.113.gtf"
bed="/home/future/GRCh38_reference/Homo_sapiens.GRCh38.113.bed"

# 创建分析目录
mkdir -p 01_fastqc_raw \
	         02_fastp_clean \
		          03_fastqc_clean \
			           04_star_align \
				            05_markdup \
					             06_rseqc \
						              07_qualimap \
							               08_rsem \
								                09_multiqc_report \
										         10_summary \
											          logs

# 输出识别到的样本ID
echo "检测到以下样本：" | tee -a logs/batch.log
for sample in "${samples[@]}"
do
	  echo "$sample" | tee -a logs/batch.log
  done
  echo "-------------------------------------------------------------" | tee -a logs/batch.log

  # 批量分析循环
  for sample in ${samples[@]}
  do
	    echo "========== Start processing sample: $sample ==========" | tee -a logs/${sample}.log

	      sample_start_time=$(date +%s)

	        echo "[FastQC] Running raw data QC..." | tee -a logs/${sample}.log
		  fastqc -t $threads \
			           -o 01_fastqc_raw \
				            00_raw_data/${sample}_1.fastq.gz \
					             00_raw_data/${sample}_2.fastq.gz 2>&1 | tee -a logs/${sample}.log
		    if [ $? -ne 0 ]; then
			        echo "ERROR: $sample FastQC 原始质控失败，退出脚本。" | tee -a logs/${sample}.log
				    exit 1
				      fi

				        echo "[Fastp] Cleaning raw reads..." | tee -a logs/${sample}.log
					  fastp -i 00_raw_data/${sample}_1.fastq.gz \
						          -I 00_raw_data/${sample}_2.fastq.gz \
							          -o 02_fastp_clean/${sample}_1.clean.fastq.gz \
								          -O 02_fastp_clean/${sample}_2.clean.fastq.gz \
									          --detect_adapter_for_pe \
										          --trim_front1 15 \
											          --trim_front2 15 \
												          --cut_tail \
													          --cut_tail_mean_quality 20 \
														          --length_required 50 \
															          --thread $threads \
																          --json 02_fastp_clean/${sample}.json \
																	          --html 02_fastp_clean/${sample}.html 2>&1 | tee -a logs/${sample}.log
					    if [ $? -ne 0 ] || [ ! -s 02_fastp_clean/${sample}_1.clean.fastq.gz ] || [ ! -s 02_fastp_clean/${sample}_2.clean.fastq.gz ]; then
						        echo "ERROR: $sample Fastp 清洗失败或输出文件缺失，退出脚本。" | tee -a logs/${sample}.log
							    exit 1
							      fi

							        echo "[FastQC] QC after cleaning..." | tee -a logs/${sample}.log
								  fastqc -t $threads \
									           -o 03_fastqc_clean \
										            02_fastp_clean/${sample}_1.clean.fastq.gz \
											             02_fastp_clean/${sample}_2.clean.fastq.gz 2>&1 | tee -a logs/${sample}.log
								    if [ $? -ne 0 ]; then
									        echo "ERROR: $sample 清洗后 FastQC 失败，退出脚本。" | tee -a logs/${sample}.log
										    exit 1
										      fi

										        echo "[STAR] Aligning reads..." | tee -a logs/${sample}.log
											  STAR --runThreadN $threads \
												         --twopassMode Basic \
													        --genomeDir $star_index \
														       --readFilesIn 02_fastp_clean/${sample}_1.clean.fastq.gz 02_fastp_clean/${sample}_2.clean.fastq.gz \
														              --readFilesCommand zcat \
															             --outFileNamePrefix 04_star_align/${sample}_ \
																            --outSAMtype BAM SortedByCoordinate \
																	           --quantMode TranscriptomeSAM \
																		          --outSAMattrRGline ID:${sample} SM:${sample} PL:ILLUMINA \
																			         --outSAMunmapped None 2>&1 | tee -a logs/${sample}.log
											    if [ $? -ne 0 ] || [ ! -s 04_star_align/${sample}_Aligned.sortedByCoord.out.bam ]; then
												        echo "ERROR: $sample STAR 比对失败或 BAM 文件缺失，退出脚本。" | tee -a logs/${sample}.log
													    exit 1
													      fi

													        echo "[samtools] Generating flagstat..." | tee -a logs/${sample}.log
														  samtools flagstat 04_star_align/${sample}_Aligned.sortedByCoord.out.bam \
															      > 04_star_align/${sample}_flagstat.txt 2>&1 | tee -a logs/${sample}.log
														    if [ $? -ne 0 ] || [ ! -s 04_star_align/${sample}_flagstat.txt ]; then
															        echo "ERROR: $sample samtools flagstat 失败或文件缺失，退出脚本。" | tee -a logs/${sample}.log
																    exit 1
																      fi

																        echo "[Picard] Marking duplicates..." | tee -a logs/${sample}.log
																	  picard MarkDuplicates \
																		           I=04_star_align/${sample}_Aligned.sortedByCoord.out.bam \
																			            O=05_markdup/${sample}_dedup.bam \
																				             M=05_markdup/${sample}_dup_metrics.txt \
																					              REMOVE_DUPLICATES=false \
																						               ASSUME_SORTED=true 2>&1 | tee -a logs/${sample}.log
																	    if [ $? -ne 0 ] || [ ! -s 05_markdup/${sample}_dedup.bam ]; then
																		        echo "ERROR: $sample Picard MarkDuplicates 失败或 BAM 缺失，退出脚本。" | tee -a logs/${sample}.log
																			    exit 1
																			      fi

																			        echo "[RSeQC] Infer experiment strand..." | tee -a logs/${sample}.log
																				  infer_experiment.py -r $bed \
																					                        -i 04_star_align/${sample}_Aligned.sortedByCoord.out.bam \
																								                      > 06_rseqc/${sample}_strand.txt 2>&1 | tee -a logs/${sample}.log
																				    if [ $? -ne 0 ] || [ ! -s 06_rseqc/${sample}_strand.txt ]; then
																					        echo "ERROR: $sample RSeQC infer_experiment 失败或结果缺失，退出脚本。" | tee -a logs/${sample}.log
																						    exit 1
																						      fi

																						        echo "[Qualimap] RNA-Seq QC report..." | tee -a logs/${sample}.log
																							  qualimap rnaseq --java-mem-size=16000M \
																								                    -bam 04_star_align/${sample}_Aligned.sortedByCoord.out.bam \
																										                      -gtf $gtf \
																												                        -outdir 07_qualimap/${sample} \
																															                  -outformat PDF:HTML \
																																	                    -pe \
																																			                      -p non-strand-specific 2>&1 | tee -a logs/${sample}.log

																							    echo "[RSEM] Calculating expression levels..." | tee -a logs/${sample}.log
																							      rsem-calculate-expression --alignments \
																								                                  --paired-end \
																												                              -p $threads \
																															                                  --strand-specific \
																																			                              --forward-prob 0 \
																																						                                  04_star_align/${sample}_Aligned.toTranscriptome.out.bam \
																																										                              $rsem_ref \
																																													                                  08_rsem/${sample} 2>&1 | tee -a logs/${sample}.log
																							        if [ $? -ne 0 ] || [ ! -s 08_rsem/${sample}.genes.results ]; then
																									    echo "ERROR: $sample RSEM 定量失败或结果缺失，退出脚本。" | tee -a logs/${sample}.log
																									        exit 1
																										  fi

																										    sample_end_time=$(date +%s)
																										      sample_elapsed=$((sample_end_time - sample_start_time))

																										        echo "样本 $sample 总耗时: ${sample_elapsed} 秒" | tee -a logs/sample_times.log
																											  echo "========== Sample $sample processing completed ==========" | tee -a logs/${sample}.log
																											    echo "样本 $sample 总耗时: ${sample_elapsed} 秒" | tee -a logs/${sample}.log
																											      echo "-------------------------------------------------------------" | tee -a logs/${sample}.log

																										      done

																										      echo "========== Generating MultiQC summary report ==========" | tee -a logs/batch.log
																										      multiqc ./ -o 09_multiqc_report 2>&1 | tee -a logs/batch.log
																										      if [ $? -ne 0 ]; then
																											        echo "ERROR: MultiQC 汇总报告失败，退出脚本。" | tee -a logs/batch.log
																												  exit 1
																										      fi

																										      echo "========== Extract mapping rate summary ==========" | tee -a logs/batch.log
																										      for i in 04_star_align/*Log.final.out
																										      do
																											        sample=$(basename $i | cut -d "_" -f1)
																												  rate=$(grep "Uniquely mapped reads %" $i | awk '{print $NF}')
																												    echo -e "${sample}\t${rate}"
																											    done > 10_summary/mapping_summary.txt

																											    echo "========== Extract duplication rate summary ==========" | tee -a logs/batch.log
																											    for i in 05_markdup/*dup_metrics.txt
																											    do
																												      sample=$(basename $i | cut -d "_" -f1)
																												        dup=$(grep "Percent Duplication" -A1 $i | tail -n1 | awk '{print $9}')
																													  echo -e "${sample}\t${dup}"
																												  done > 10_summary/duplication_summary.txt

																												  echo "========== All batch processing finished ==========" | tee -a logs/batch.log

																												  # 生成项目综合 summary 报告
																												  echo "========== Generating Final Analysis Summary ==========" | tee -a logs/batch.log

																												  summary_file="README.summary"
																												  echo "================= RNA-seq 批量分析执行记录 =================" > $summary_file
																												  echo "分析日期       : $(date)" >> $summary_file
																												  echo "执行机器       : $(hostname)" >> $summary_file
																												  echo "" >> $summary_file

																												  echo "分析样本数量   : ${#samples[@]}" >> $summary_file
																												  echo "样本ID列表     :" >> $summary_file
																												  for sample in "${samples[@]}"
																												  do
																													    echo "  - $sample" >> $summary_file
																												    done
																												    echo "" >> $summary_file

																												    echo "软件版本：" >> $summary_file
																												    echo "  fastp      : $(fastp --version | head -n1)" >> $summary_file
																												    echo "  STAR       : $(STAR --version | head -n1)" >> $summary_file
																												    echo "  samtools   : $(samtools --version | head -n1 | awk '{print $2}')" >> $summary_file
																												    echo "  RSEM       : $(rsem-calculate-expression --version 2>&1 | head -n1)" >> $summary_file
																												    echo "  MultiQC    : $(multiqc --version)" >> $summary_file
																												    echo "" >> $summary_file

																												    echo "样本运行耗时（单位：秒）：" >> $summary_file
																												    cat logs/sample_times.log >> $summary_file
																												    echo "" >> $summary_file

																												    echo "Mapping Rate 汇总：" >> $summary_file
																												    cat 10_summary/mapping_summary.txt >> $summary_file
																												    echo "" >> $summary_file

																												    echo "Duplication Rate 汇总：" >> $summary_file
																												    cat 10_summary/duplication_summary.txt >> $summary_file
																												    echo "" >> $summary_file

																												    echo "输出目录结构（层级2）：" >> $summary_file
																												    tree -L 2 . >> $summary_file
																												    echo "" >> $summary_file

																												    echo "全部日志文件路径：" >> $summary_file
																												    echo "  logs/ 目录内为单样本分析日志和批量分析日志" >> $summary_file
																												    echo "" >> $summary_file

																												    echo "================= 执行记录生成完毕 =================" | tee -a logs/batch.log
																												    echo "最终综合记录保存在: $summary_file" | tee -a logs/batch.log
