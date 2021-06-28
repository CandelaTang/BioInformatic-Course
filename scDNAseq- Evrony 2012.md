## Single-cell DNA seq
### Reproduce from Gilad D. Evrony et al., Cell, 2012
#### Codes
1. Bamcoverge
```json
#For homework 3
#bamcoverage.py

#get bins from reference sample: 1465-cortex_50k-presort_MDA_5_WGS
#consolidate bins to get bins with equal number of reads

cmd = '''
bamCoverage --bam {bamfile}.sorted.bam -o {bamfile}.{binsize}.bdg \
    --binSize {binsize} \
    --outFileFormat bedgraph && \
grep chr {bamfile}.{binsize}.bdg | grep -v "chrM" {bamfile}.{binsize}.bdg > {bamfile}.{binsize}.c24.bdg && \
ln -s {bamfile}.{binsize}.c24.bdg
'''

#拆分为两步进行
#先做下第一步bamcoverage

import os

cmd = '''
bamCoverage --bam {bamfile}.sorted.bam -o {bamfile}.{binsize}.tky.bdg \
    --binSize {binsize} \
    --outFileFormat bedgraph
'''


#reference: SRR536742 1465-cortex_50k-presort_MDA_1_WGS
#sample: SRR536764 866-cortex_1-neuron_MDA_2_WGS
#use binned result for SRR536742 by CSY, only need to process SRR536764 here
s = [ 'SRR536764' ]
binsize = 5000

p2data = '/data/user_05/projects/gilad_cell_2012/data'  #enter your sra data directory

shtxt = []
for name in s:
	p2bam = os.path.join(os.path.join(p2data,name),name)
	cmdstr = cmd.format(bamfile=p2bam,binsize=binsize)
	print(cmdstr)
	shtxt.append(cmdstr)

with open( 'tkybamcov.sh' , 'w' ) as fh:
    fh.write( '#!/usr/bin/bash\n' )
    fh.write( '#SBATCH -J bamcovtky\n') #queue in cluster
    fh.write( '#SBATCH -p CN_BIOT\n')
    fh.write( '#SBATCH --nodes=1\n')
    fh.write( '#SBATCH --ntasks=4\n')
    fh.write( '#SBATCH --output=%j.out\n')
    fh.write( '#SBATCH --error=%j.err\n')
    fh.write( '\n'.join(shtxt) )
    

#enter $sbatch tkybamcov.sh to run


#再做第二步提取24个chr的数据

import os
cmd = '''
grep chr {bamfile}.{binsize}.tky.bdg | grep -v "chrM" {bamfile}.{binsize}.tky.bdg > {bamfile}.{binsize}.c24.tky.bdg && \
ln -s {bamfile}.{binsize}.c24.tky.bdg
'''

s = [ 'SRR536764' ]
binsize = 5000

p2data = '/data/user_05/projects/gilad_cell_2012/data'  #enter your sra data directory

shtxt = []
for name in s:
	p2bam = os.path.join(os.path.join(p2data,name),name)
	cmdstr = cmd.format(bamfile=p2bam,binsize=binsize)
	print(cmdstr)
	shtxt.append(cmdstr)

with open( 'tkybamcov2.sh' , 'w' ) as fh:
    fh.write( '#!/usr/bin/bash\n' )
    fh.write( '#SBATCH -J bamcovtky2\n') #queue in cluster
    fh.write( '#SBATCH -p CN_BIOT\n')
    fh.write( '#SBATCH --nodes=1\n')
    fh.write( '#SBATCH --ntasks=4\n')
    fh.write( '#SBATCH --output=%j.out\n')
    fh.write( '#SBATCH --error=%j.err\n')
    fh.write( '\n'.join(shtxt) )
    
#$sbatch tkybamcov2.sh to run

#go to getbins.m
```
2. get bins for genome
```json
%getbins.m
% get bins for genome
a = importdata('SRR536742.5000.c24.bdg'); %enter reference bdg file
dat = a.data;
chr_names = a.textdata(1:end,1);

%start, end, counts
%{
get total counts
get counts per bin = total counts divided by 6,000
for each chromosome
    get #bins per chromosome, take ceiling
    get coordinates for each bin
%}

%get total counts
tc = sum(dat(:,3));
nbin=6e3; %create 6,000 bins
cpb = tc/nbin;
%get index for the first segment of each chromosome
i_chr1seg = find([dat(1:end,1)]<[inf;dat(1:end-1,1)]);
chr_names(i_chr1seg)
dat(i_chr1seg,:)
i_chr1seg = [i_chr1seg;size(dat,1)+1];

%for each chromosome
binbychr = {};
for i = 1:length(i_chr1seg)-1
    %get #bins per chromosome, take ceiling
    i_chr_start = i_chr1seg(i);
    i_chr_end = i_chr1seg(i+1)-1;
    dat_chr = dat(i_chr_start:i_chr_end,:);
    name_chr = chr_names(i_chr_start:i_chr_end,:);
    tc_chr = sum(dat_chr(:,3));
    nbin_chr = ceil(tc_chr/cpb);
    
    %get coordinates for each bin
    x = cumsum(dat_chr(:,3))/tc_chr;
    y = linspace(0,1,nbin_chr)';
    iy = zeros(size(y));iy(1)=1;iy(end) = length(x)+1;
    for j = 2:length(y)-1
        iy(j) = find(x>y(j),1);
    end
    
    %format start, end, counts for each bin
    bin_start = dat_chr(iy(1:end-1),1);
    bin_end = [bin_start(2:end);dat_chr(iy(end)-1,2)];
    bin_count = zeros(size(bin_start));
    for j = 1:length(iy)-1
        bin_count(j) = sum(dat_chr(iy(j):iy(j+1)-1,3));
    end
    
    %save results
    binbychr{i,1}=name_chr(iy(1:end-1));
    binbychr{i,2}=[bin_start,bin_end,bin_count];
end

%write to bed graph file
fid = fopen('chr_count_bin.bdg','w');
for i = 1:size(binbychr,1)
    chrtext = binbychr{i,1};
    chrdata = binbychr{i,2};
    for j = 1:size(chrtext,1)
        fprintf(fid,'%s\t%d\t%d\t%d\n',chrtext{j},chrdata(j,1),chrdata(j,2),chrdata(j,3));
    end
end
fclose(fid);
```
3. multibamsummary
```json
#multibamsummary.py
#after getbins.m
#count B_Ref, B_Sample

#get coverage for bins
import os
cmd = '''
multiBamSummary BED-file --BED {bedfile} --bamfiles {bamfiles} -o \
    {outputnpz}.npz --labels {samplenames} --outRawCounts {outputtab}.tab
    '''
#enter your own directory
p2data = '/data/user_05/projects/gilad_cell_2012/data'
sratable = '/data/user_05/projects/gilad_cell_2012/results/tangky/SraRunTable.tky.csv'

#enter your reference and sample sra
sra_list = '''
SRR536742
SRR536742
SRR536743
SRR536750
SRR536751
SRR536757
SRR536758
SRR536763
SRR536764
'''
sra_list = sra_list.strip().split('\n')

srr2name = dict()
with open(os.path.join(p2data,sratable),'r') as fh:
    for line in fh:
        srracc, samplename = line.strip().split(',')[:2]
        srr2name[srracc]=samplename
        

p2bamlist_str = []
labels_str = []
for srracc in sra_list:
    p2bam = os.path.join(os.path.join(p2data,srracc),srracc) + '.sorted.bam'
    p2bamlist_str.append(p2bam)
    labels_str.append(srr2name[srracc])
    
p2bamlist_str =' '.join(p2bamlist_str)
labels_str =' '.join(labels_str)

outputnpz = '/data/user_05/projects/gilad_cell_2012/results/tangky/counts_by_bins0528' #enter your output directory
outputtab = outputnpz
bedfile = '/data/user_05/projects/gilad_cell_2012/results/tangky/chr_count_bin.bdg' #enter your own bdg directory

shtxt = []
cmdstr = cmd.format(bedfile=bedfile, bamfiles=p2bamlist_str,
                    outputnpz=outputnpz, samplenames=labels_str, outputtab=outputtab)
print(cmdstr)
shtxt.append(cmdstr)
with open( 'tkymtbam0528.sh' , 'w' ) as fh:
    fh.write( '#!/usr/bin/bash\n' )
    fh.write( '#SBATCH -J mtbamtky\n') #queue in cluster
    fh.write( '#SBATCH -p CN_BIOT\n')
    fh.write( '#SBATCH --nodes=1\n')
    fh.write( '#SBATCH --ntasks=4\n')
    fh.write( '#SBATCH --output=%j.out\n')
    fh.write( '#SBATCH --error=%j.err\n')
    fh.write( '\n'.join(shtxt) )
    
#enter $sbatch tkymtbam.sh to run
#go to countbins.m
```
4. count and process bins
```json
%countbins.m
%count and process bins

%import bin coverage
a = importdata('counts_by_bins0528.tab'); %enter your tab file
dat = double(a.data);
sample_names = a.textdata(1,2:end);
for i =1:length(sample_names)
    sample_names{i} = strrep(sample_names(i),"'","");
end
chr_names = a.textdata(2:end,1);

%get index for the first segment of each chromosome
i_chr1seg = find([dat(1:end,1)]<[inf;dat(1:end-1,1)]);
chr_names(i_chr1seg)
dat(i_chr1seg,:)
i_chr1seg = [i_chr1seg;size(dat,1)+1];

%%calculate
%calculate CN_n = (B_sample,n/T_sample)/(B_ref/T_ref)
dat_cn = zeros(size(dat));
T = sum(dat,1);
i_ref = 3; %identify your ref
for i = 3:size(dat,2)
    dat_cn(:,i)=(dat(:,i)./T(i))./(dat(:,i_ref)./(T(i_ref)));
end

%calculate median of B_sample,n for copynumber > 0.5 (median of all
%samples)
B_med = zeros(size(dat,1),1);
for j = 1:size(dat,1)
    col = find(dat_cn(j,:)>0.5);
    B_med(j) = median(dat(j,col));
end

%normalize B_sample,n by median
cn_log = zeros(size(dat));
cn_log(:,1:2) = dat(:,1:2);
for j = 1:size(dat,1)
    cn_log(j,3:end) = log2(dat(j,3:end)./B_med(j));
end


i_sample = [8 9 10 11]; %enter your sample


%%PLOT
figure;
for i = 1:length(i_sample)
    subplot(length(i_sample),1,i)
    plot(cn_log(:,i_sample(i)),'.');
    hold on
    for j = 2:length(i_chr1seg)-1
        line([i_chr1seg(j),i_chr1seg(j)],[4,-4],...
            'linewidth',0.5,...
            'color',[1 0 0])
    end
    plot([1,size(cn_log,1)],[1,1],'--m')
    plot([1,size(cn_log,1)],[0,0],'--k')
    plot([1,size(cn_log,1)],[-1,-1],'--m')
    
    ylabel(sample_names(i_sample(i)),'Interpreter','none');
    xtck = (i_chr1seg(1:end-1) + i_chr1seg(2:end))/2;
    xtickslbl = cellstr(num2str([1:24]'));
    xtickslbl{23} = 'X';
    xtickslbl{24} = 'Y';
    set(gca,'xtick',xtck)
    set(gca,'xticklabel',xtickslbl)
    set(gca,'ylim',[-4,4])
end
```
