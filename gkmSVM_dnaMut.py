# f = open( 'test.txt' , 'r' )
# seq_line = f.readline()
# dna_seq = ["A", "T", "G", "C"]


pwd = '/Users/suryachhetri/Dropbox/for_tools'
dna_seq = ["A", "T", "G", "C"]

#in silico mutagenesis of full sequence
#mutate each postion with A/T/G/C for all possible combination scores
#acrss the sequence
dna_seq = set("ATGC")
with open('test.txt') as line:
	for idx, value in enumerate(line):
		seq_value = value.rstrip()
		#print(seq_value)
		if not seq_value.lstrip().startswith('>'):
			print(seq_value)
			seq_count=len(seq_value)
			#print(seq_count)
			for i in range(seq_count):
				base = seq_value[i].upper()
                print("\n")
				print("base position: {} i.e., {} > .".format(i, base))
                #print reference base sequence
                print("{} \t > base {}:{}|{} --> \t {}".format(i, i, base, base, seq_value))

				seq = set(base)
				unique_seq = dna_seq - seq
				for char in unique_seq:
					position = i
					new_character = char
					string = seq_value[:position] + new_character + seq_value[position+1:]
					#print("old {} _ {}".format(idx, seq_value))
                    print("{} \t > base {}:{}>{} --> \t {}".format(position, position, base, char, string))


#check the impact of fixed varaint "T" across the given sequence
#mutate only single "fixed" base across the sequence, i.e., T in this case
mutate = "C"
dna_seq = set("ATGC")
with open('test.txt') as line:
    for idx, value in enumerate(line):
        seq_value = value.rstrip()
        #print(seq_value)
        if not seq_value.lstrip().startswith('>'):
            print(seq_value)
            seq_count=len(seq_value)
            #print(seq_count)
            # for i in range(seq_count):
            #     base = seq_value[i].upper()
            #     print("\n")
            #     print("base position: {} i.e., {} > .".format(i, base))
            #     print("{} \t > base {}:{}|{} --> \t {}".format(i, i, base, base, seq_value))


            #seq = set(base)
            #unique_seq = dna_seq - seq
            print("\nI. \t > base ref --> \t {}\n".format(seq_value))
            for idx, value in enumerate(list(seq_value)):
                position = idx
                new_character = mutate
                string = seq_value[:position] + new_character + seq_value[position+1:]
                
                #print("old {} _ {}".format(idx, seq_value))
                print("{} \t > base {}:{}>{} --> \t {}".format(position, position, value, new_character, string))



pwd = '/Users/suryachhetri/Dropbox/for_tools'



#1) Input SNPs bed format: chr2:16:T>A (#eQTLs/GWAS SNPs)
# chr   start   end    ref  alt
# chr 2   15    16      T   A 

#2) Flank sequence 10bp each side (Bedtools)
# Chrom   start   end    ref  alt
# chr 2   6      26      T   A


#3)Run getfasta Bedtools #Retain SNP-name chr2:16:T>A
#> chr2:6-26
#ACCCGTGGGT|T|AGGGCCCTAT


# '''BEDTOOLS GETFASTA Custom name: create 4th cols'''

# # -name Using the BED “name” column as a FASTA header.
# # Using the -name option, one can set the FASTA header 
# #for each extracted sequence to be the “name” columns 
# #from the BED feature.

# # $ cat test.fa
# # >chr1
# # AAAAAAAACCCCCCCCCCCCCGCTACTGGGGGGGGGGGGGGGGGG

# #Create 4th column
## $7 is snps position
# #awk -v OFS="\t" '{$4=$1":"$2"-"$3":"$7"-"$4"-"$5":"$6}1' test.bed
# # $ cat test.bed
# # chr1 5 10 chr1:6-26:16-T-A:TFname

# # $ bedtools getfasta -fi test.fa -bed test.bed -name
# # >chr1:5-10:T:A_GENE(TF-NAME)
# # AAACC

# #It seems name comes from bed. 
# #If so you can create a fifth column with desired string.

# # awk -v OFS="\t" '{$5=$4":"$1":"$2"-"$3}1' test.bed | bedtools getfasta -fi test.fa -bed - -name
# # >Name:chr1:5-10
# # AAACC


#4) Read SNPs list containing fasta file

#def read_file():
#Read 2 lines at a time: zip(f,f)
#Read 3 lines at a time: for l1,l2,l3 in zip(f,f,f):
# with open('test_snpfasta.txt') as f:
#     for line1,line2 in zip(f,f):
#         print(line1,line2)
#         print("Fasta change ...")

#import itertools
#with open('test_snpfasta.txt') as f:
#with open('test_snpfasta_new.txt'), open('ref.fa', 'w'), open('alt.fa', 'w') as f, ref_file, alt_file:

#5) Write Fasta containing file with base edits


###################################################

# FINALIZED SCRIPT STARTS HERE:

###################################################

import itertools

# Single point mutation and not across the sequence
with open('ref.fa', 'w') as ref_file:
    with open('alt.fa', 'w') as alt_file:
        #snps list containing fasta file derived above getfasta
        with open('test_snpfasta_new.txt') as f:

            for line1,line2 in itertools.zip_longest(*[f]*2):
                print("\n\nFasta SEQ ...")
                print(line1,line2)

                #remove \n characters and find snp pos
                #locus_pos = line1.rstrip("\n").split(":")[1].split("-")
                chr_pos = line1.rstrip("\n").split(":")[0]
                print("Chr: {}".format(chr_pos))
                snp_id = line1.rstrip("\n").split(":")[2]
                print("SNP Id: {}".format(snp_id))
                locus_pos = line1.rstrip("\n").split(":")[1].split("-")
                print("Locus postion: {}".format(locus_pos))
                snp_edit_pos = int((int(locus_pos[1]) - int(locus_pos[0]))/2)
                print("SNP edit postion: {}th".format(snp_edit_pos))
                
                base_mutate = line1.rstrip("\n").split(":")[2].split("-")[2:]
                print("Base mutate: {}".format(base_mutate))
                tf_name = line1.rstrip("\n").split(":")[3]
                print("Annotation: {}\n".format(tf_name))

                #Load fasta sequence
                ori_fasta = line2.rstrip("\n")

                #SNP mutate
                mutate_chars = line1.rstrip("\n").split(":")[2].split("-")[2:]
                position = snp_edit_pos
                ref_base = mutate_chars[0]
                mut_base = mutate_chars[1]
                mut_fasta = ori_fasta[:position] + mut_base + ori_fasta[position+1:]
                print("Ori Sequence: {}".format(ori_fasta))
                print("Mut Sequence: {}".format(mut_fasta))

                #Write output files
                ref_file.write(line1.rstrip("\n")+"_"+"ref"+"\n")
                ref_file.write(ori_fasta+"\n")
                alt_file.write(line1.rstrip("\n")+"_"+"alt"+"\n")
                alt_file.write(mut_fasta+"\n")

                # ref_file.write(line1)
                # ref_file.write(ori_fasta+"\n")
                # alt_file.write(line1)
                # alt_file.write(mut_fasta+"\n")

        #print("old {} _ {}".format(idx, seq_value))
        #print("{} \t > base {}:{}>{} --> \t {}".format(position, position, value, char, string))




#in silico mutagenesis of full sequence
#mutate each postion with A/T/G/C for all possible combination scores
#acrss the sequence
dna_seq = set("ATGC")
# with open('ref.fa', 'w') as ref_file:
with open('base_mutate_saturated_final.fa', 'w') as mut_file:
    #snps list containing fasta file derived above getfasta
    with open('test_snpfasta_new.txt') as f:
        for line1,line2 in itertools.zip_longest(*[f]*2):
            print("\n\nFasta SEQ ...")
            print(line1,line2)

            #remove \n characters and find snp pos
            #locus_pos = line1.rstrip("\n").split(":")[1].split("-")
            chr_pos = line1.rstrip("\n").split(":")[0]
            print("Chr: {}".format(chr_pos))
            snp_id = line1.rstrip("\n").split(":")[2]
            print("SNP Id: {}".format(snp_id))
            locus_pos = line1.rstrip("\n").split(":")[1].split("-")
            print("Locus postion: {}".format(locus_pos))
            snp_edit_pos = int((int(locus_pos[1]) - int(locus_pos[0]))/2)
            print("SNP edit postion: {}th".format(snp_edit_pos))
            
            base_mutate = line1.rstrip("\n").split(":")[2].split("-")[2:]
            print("Base mutate: {}".format(base_mutate))
            tf_name = line1.rstrip("\n").split(":")[3]
            print("Annotation: {}\n".format(tf_name))


            seq_value = line2.rstrip("\n")
            print(seq_value)
            seq_count=len(seq_value)
            #print(seq_count)
            for i in range(seq_count):
                base = seq_value[i].upper()
                print("\n")
                print("base position: {} i.e., {} > .".format(i, base))
                #print reference base sequence
                print("{} \t > base {}:{}|{} --> \t {}".format(i, i, base, base, seq_value))

                #Write base reference output files
                #mut_file.write(line1)
                mut_file.write(line1.rstrip("\n")+"_"+str(i)+base+"\n")
                mut_file.write(seq_value+"\n")

                seq = set(base)
                unique_seq = dna_seq - seq
                for char in unique_seq:
                    position = i
                    new_character = char
                    string = seq_value[:position] + new_character + seq_value[position+1:]
                    #print("old {} _ {}".format(idx, seq_value))
                    print("{} \t > base {}:{}>{} --> \t {}".format(position, position, base, char, string))

                    #print("Ori Sequence: {}".format(ori_fasta))
                    #print("Mut Sequence: {}".format(mut_fasta))

                    #Write output files
                    #mut_file.write(line1)
                    mut_file.write(line1.rstrip("\n")+"_"+str(i)+str(char)+"\n")
                    mut_file.write(string+"\n")
 


#check the impact of fixed varaint "T" across the given sequence
#mutate only single "fixed" base across the sequence, i.e., T in this case
#mutate = "C"
dna_seq = set("ATGC")
with open('base_mutate_fixedSNP_final.fa', 'w') as mut_file:
    #snps list containing fasta file derived above getfasta
    with open('test_snpfasta_new.txt') as f:
        for line1,line2 in itertools.zip_longest(*[f]*2):
            print("\n\nFasta SEQ ...")
            print(line1,line2)

            #remove \n characters and find snp pos
            #locus_pos = line1.rstrip("\n").split(":")[1].split("-")
            chr_pos = line1.rstrip("\n").split(":")[0]
            print("Chr: {}".format(chr_pos))
            snp_id = line1.rstrip("\n").split(":")[2]
            print("SNP Id: {}".format(snp_id))
            locus_pos = line1.rstrip("\n").split(":")[1].split("-")
            print("Locus postion: {}".format(locus_pos))
            snp_edit_pos = int((int(locus_pos[1]) - int(locus_pos[0]))/2)
            print("SNP edit postion: {}th".format(snp_edit_pos))
            
            base_mutate = line1.rstrip("\n").split(":")[2].split("-")[2:]
            print("Base mutate: {}".format(base_mutate))
            tf_name = line1.rstrip("\n").split(":")[3]
            print("Annotation: {}\n".format(tf_name))


            seq_value = line2.rstrip("\n")
            print(seq_value)
            seq_count=len(seq_value)
            print("\nI. \t > base ref --> \t {}\n".format(seq_value))

            #Write base reference output files
            #mut_file.write(line1)
            mut_file.write(line1.rstrip("\n")+"_"+"alt"+"\n")
            mut_file.write(seq_value+"\n")

            #SNP mutate
            mutate_chars = line1.rstrip("\n").split(":")[2].split("-")[2:]
            ref_base = mutate_chars[0]
            mut_base = mutate_chars[1]


            #print("\nI. \t > base ref --> \t {}\n".format(seq_value))
            for idx, value in enumerate(list(seq_value)):
                position = idx
                new_character = mut_base
                string = seq_value[:position] + new_character + seq_value[position+1:]
                
                #print("old {} _ {}".format(idx, seq_value))
                print("{} \t > base {}:{}>{} --> \t {}".format(position, position, value, new_character, string))

                #Write output files
                #mut_file.write(line1)
                mut_file.write(line1.rstrip("\n")+"_"+str(idx)+mut_base+"\n")
                mut_file.write(string+"\n")




# Single point mutation plus sliding window and not across the sequence
with open('ref_slide_basemut.fa', 'w') as ref_file:
    with open('alt_slide_basemut.fa', 'w') as alt_file:
        with open('test_snpfasta_new.txt') as f:

            for line1,line2 in itertools.zip_longest(*[f]*2):
                print("\n\nFasta SEQ ...")
                print(line1,line2)

                #remove \n characters and find snp pos
                #locus_pos = line1.rstrip("\n").split(":")[1].split("-")
                chr_pos = line1.rstrip("\n").split(":")[0]
                print("Chr: {}".format(chr_pos))
                snp_id = line1.rstrip("\n").split(":")[2]
                print("SNP Id: {}".format(snp_id))
                locus_pos = line1.rstrip("\n").split(":")[1].split("-")
                print("Locus postion: {}".format(locus_pos))
                snp_edit_pos = int((int(locus_pos[1]) - int(locus_pos[0]))/2)
                print("SNP edit postion: {}th".format(snp_edit_pos))
                
                base_mutate = line1.rstrip("\n").split(":")[2].split("-")[2:]
                print("Base mutate: {}".format(base_mutate))
                tf_name = line1.rstrip("\n").split(":")[3]
                print("Annotation: {}\n".format(tf_name))

                #Load fasta sequence
                ori_fasta = line2.rstrip("\n")

                #SNP mutate
                mutate_chars = line1.rstrip("\n").split(":")[2].split("-")[2:]
                position = snp_edit_pos
                ref_base = mutate_chars[0]
                mut_base = mutate_chars[1]
                mut_fasta = ori_fasta[:position] + mut_base + ori_fasta[position+1:]
                slide_ref_fasta = ori_fasta[:position+1]
                slide_alt_fasta = ori_fasta[:position] + mut_base
                print("\nOri Sequence: {} \t \t Mut Sequence: {}\n".format(ori_fasta, mut_fasta))
                #print("Mut Sequence: {}".format(mut_fasta))
                #print("Ori Sequence: {}".format(slide_ref_fasta))
                #print("Mut Sequence: {}".format(slide_alt_fasta))
                
                #Write output original fasta files for reference
                ref_file.write(line1.rstrip("\n")+"_"+"ref"+"\n")
                ref_file.write(ori_fasta+"\n")
                alt_file.write(line1.rstrip("\n")+"_"+"alt"+"\n")
                alt_file.write(mut_fasta+"\n")

                # ref_file.write(line1)
                # ref_file.write(ori_fasta+"\n")
                # alt_file.write(line1)
                # alt_file.write(mut_fasta+"\n")

                
                for pos in range(position+1):
                    slide_ref_fasta = ori_fasta[pos:position+1+pos]
                    slide_alt_fasta = mut_fasta[pos:position+1+pos]
                    print("{} \t \t {}".format(slide_ref_fasta, slide_alt_fasta))

                    #Write output files
                    ref_file.write(line1.rstrip("\n")+"_"+str(pos)+"\n")
                    ref_file.write(slide_ref_fasta+"\n")
                    alt_file.write(line1.rstrip("\n")+"_"+str(pos)+"\n")
                    alt_file.write(slide_alt_fasta+"\n")

                    # ref_file.write(line1)
                    # ref_file.write(slide_ref_fasta+"\n")
                    # alt_file.write(line1)
                    # alt_file.write(slide_alt_fasta+"\n")



# In [59]: dna_seq = set("ATGC")
#     ...: with open('test1.txt') as line:
#     ...: ^Ifor idx, value in enumerate(line):
#     ...: ^I^Iseq_value = value.rstrip()
#     ...: ^I^I#print(seq_value)
#     ...: ^I^Iif not seq_value.lstrip().startswith('>'):
#     ...: ^I^I^Iprint(seq_value)
#     ...: ^I^I^Iseq_count=len(seq_value)
#     ...: ^I^I^I#print(seq_count)
#     ...: ^I^I^Ifor i in range(seq_count):
#     ...: ^I^I^I^Ibase = seq_value[i].upper()
#     ...: ^I^I^I^I#print(i)
#     ...:
#     ...: ^I^I^I^Iseq = set(base)
#     ...: ^I^I^I^Iunique_seq = dna_seq - seq
#     ...: ^I^I^I^Ifor char in unique_seq:
#     ...: ^I^I^I^I^Iposition = i
#     ...: ^I^I^I^I^Inew_character = char
#     ...: ^I^I^I^I^Istring = seq_value[:position] + new_character + seq_value[position+1:]
#     ...: ^I^I^I^I^I#print("old {} _ {}".format(idx, seq_value))
#     ...: ^I^I^I^I^Iprint("new {} _ {}".format(idx, string))