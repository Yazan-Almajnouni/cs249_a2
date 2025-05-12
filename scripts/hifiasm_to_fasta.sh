# awk '/^S/{print ">"$2;print $3}' /ibex/user/majnouym/cs249_a2/lizard_asembly/hifiasm/hifiasm_lizard.bp.p_ctg.gfa > /ibex/user/majnouym/cs249_a2/lizard_asembly/hifiasm/hifiasm_asm.fa

awk '/^S/{print ">"$2;print $3}' /ibex/user/majnouym/cs249_a2/new_assembly/filtered_reads/hifiasm/hifiasm_lizard.bp.p_utg.gfa > /ibex/user/majnouym/cs249_a2/new_assembly/filtered_reads/hifiasm/new_asm.fasta