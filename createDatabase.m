clearvars;

% split the file 'Blood.gbk' into 6561 frames under frames_folder
frames_folder = 'HMP_frames_oral/';
directory = dir(frames_folder); % names and sizes of the files
frames = char(directory.name); % File names with .txt extension 
frames = frames(frames(:,1)=='f',1:size(frames,2)); % frame file names
num_frame = length(frames); % number of frames
num_protein = 0; % number of proteins total in dataset
num_bsEnz = 0; % number of biosynthetic enzymes
%aa_counts = zeros(20,1); % counts number of enzymes for each amino acid
org_count = 0;
%scoref = [];
bcheck = 0;
threshold = 1e-5;
scoref = NaN(1e7,1);
prots = 0;


orgData = struct([]);
blastData = struct('CDS',[]);

blastSeqs = load('Tables/blast_seq_info.mat');
aminos = blastSeqs.aa;
seqList = blastSeqs.seq;

% process each frame
%for idx = 1:300
for idx = 1:num_frame
    
    % progress update
	disp([num2str(idx/num_frame*100,3) '% complete. '  num2str(num_bsEnz) ' biosynthetic enzymes found.' ]);
    %display(idx)
    
    % use genbankread function to read frame
    try
        data = genbankread([frames_folder 'frame' num2str(idx) '.txt']);
    catch
        continue;
    end
    
    % creat a new organism data for the first frame or a new organism
    if org_count == 0 || ~strcmp(orgData(org_count).name,data.Source)
        added = 0; % total number of protein
        BS_added = 0;
        totalLength = 0; % total length of the protein
        totalNum = zeros(1,20);
        org_count = org_count + 1;
        
        orgData(org_count).name = data.Source;
        orgData(org_count).BSgene.n = zeros(1,20);
        orgData(org_count).BSgene.count = [];
        orgData(org_count).BSgene.product = '';
        orgData(org_count).BSgene.l_enz = 0;
        %orgData(org_count).BSgene.CB = nan(1,20);
        %orgData(org_count).BSgene.CB_prev = nan(1,20);
        orgData(org_count).enzNum = 0;
        orgData(org_count).l_avg = 0;
        orgData(org_count).p = zeros(1,20);
        orgData(org_count).bodysite = 'Air';
        orgData(org_count).CB = inf(1,20);
        orgData(org_count).CB_prev = inf(1,20);
        
        taxon_id = findTaxon(data);
        orgData(org_count).taxon_id = taxon_id;
        
        % add one case that only two rows in SourceOrganism
        % second row only four ; and one .
        % add the strings 2:end together using + and then check
        [m,n] = size(data.SourceOrganism);
        info = string(data.SourceOrganism(m-1,:)) + ' ' + string(data.SourceOrganism(m,:));
        % five ; one .
        % four ; one .
        % special cases: oral frame 1 and oral frame 501
        loc = strfind(info,';');
        loc_end = strfind(info,'.');
        info = char(info);
        if length(loc) == 5
            orgData(org_count).kingdom = info(1:loc(1)-1);
            orgData(org_count).phylum = info(loc(1)+2:loc(2)-1);
            orgData(org_count).class = info(loc(2)+2:loc(3)-1);
            orgData(org_count).order = info(loc(3)+2:loc(4)-1);
            orgData(org_count).family = info(loc(4)+2:loc(5)-1);
            orgData(org_count).genus = info(loc(5)+2:loc_end(length(loc_end))-1);
        end
        if length(loc) == 4
            orgData(org_count).kingdom = info(1:loc(1)-1);
            orgData(org_count).phylum = info(loc(1)+2:loc(2)-1);
            orgData(org_count).class = info(loc(2)+2:loc(3)-1);
            orgData(org_count).order = info(loc(3)+2:loc(4)-1);
            orgData(org_count).family = info(loc(4)+2:loc_end(length(loc_end))-1);
            orgData(org_count).genus = '';
        end
        
        orgData(org_count).high_quality = 0;
        if strfind(data.Keywords,'HIGH')
            orgData(org_count).high_quality = 1;
        end
    end
    
    % to check CDS exist
    if isfield(data,'CDS')
        for i = 1:length(data.CDS)
            added = added + 1; % update the number of protein in current organism
            
            % use EC number to decide whether it is a biosynthetic enzyme
            [check,count] = isBiosyntheticEnzymeEC(data.CDS(i).EC_number);
            if count == 4
                disp('++found++')
            end
            % if check is 0, use gene to decide if it is a biosynthetic enzyme
            
            % Next few lines are for debugging and should be removed, all
            % this does is return the scores of the alignment to see what
            % values we can expect for the alignment score if there is a
            % match or no match
            %[check2,score] = isBiosyntheticEnzymeBlast(data.CDS(i).translation); 
            %scoref = [scoref score];            
            %if score > 1
            %    display(['ENZYME FOUND!!!!  CHECK = ' num2str(count) '   SCORE = ' num2str(score)])
            %end
            
            if ~check
                
                [check,count] = isBiosyntheticEnzymeEC(data.CDS(i).gene);
                
                if ~check % if no match with EC number and gene name, check
                 % sequence alignment
                
                if isempty(data.CDS(i).translation)
                    disp('EMPTY SEQUENCE')
                    continue;
                end
                 
                [check,count,score] = isBiosyntheticEnzymeBlastNew(data.CDS(i).translation,seqList,aminos,threshold); 
                %prots = prots+1;
                %scoref(prots) = max(score);
                
                
                if check
                    disp(['Blasted check = ' num2str(count)])
                    bcheck = bcheck + 1;
                    blastData(bcheck).CDS = data.CDS(i);
                    
                end
                
                end
            end

            seq = data.CDS(i).translation;
            seqLength = length(seq);
            totalLength = totalLength + seqLength;
            if check
                % update number of BS genes
                BS_added = BS_added + 1;
                % update count for protein
                orgData(org_count).BSgene(BS_added).count = count;
                % update product of BS gene
                orgData(org_count).BSgene(BS_added).product = data.CDS(i).product;
                % update length 
                orgData(org_count).BSgene(BS_added).l_enz = seqLength;
                orgData(org_count).BSgene(BS_added).CDS = data.CDS(i);
                
                orgData(org_count).enzNum = orgData(org_count).enzNum + 1;
                num_bsEnz = num_bsEnz + 1;
            end
            

            if check
                end
            % update n for protein
            for j = 1:20
                
                str = aminolookup('INTEGER',j);% use aminolookup and get the first letter
                enz = str(1);
                % count the number of this amino acid in this protein
                num = 0;
                
                for k = 1:seqLength
                    
                    if seq(k) == enz
                         
                         num = num + 1;
                         totalNum(j) = totalNum(j) + 1;
                         
                    end
                    
                end
                
                % save num for BS genes
                if check
                    orgData(org_count).BSgene(BS_added).n(j) = num;
                end
                
            end
            
            % update average length of all protein
            orgData(org_count).l_avg = totalLength / added;
            
            for j = 1:20
                % update percentage of cognate amino acid
                orgData(org_count).p(j) = totalNum(j) / totalLength;
            end
        end
        %num_protein = num_protein + length(data.CDS); % Track total number of protein
    else
        fprintf('No CDS present in frame %d.', idx);
        continue;
    end
end

disp(['100% complete. '  num2str(num_bsEnz) ' biosynthetic enzymes found.' ]);

% compute cognate bias and length independent cognate bias
for idx = 1:length(orgData)
    for i = 1:length(orgData(idx).BSgene)
        for j = 1:20
            orgData(idx).BSgene(i).CB(j) = NaN;
            orgData(idx).BSgene(i).CB_prev(j) = NaN;
        end
        for j = 1:length(orgData(idx).BSgene(i).count)
            aa = orgData(idx).BSgene(i).count(j);
            CB = orgData(idx).BSgene(i).n(aa)-orgData(idx).p(aa)*orgData(idx).BSgene(i).l_enz;
            CB_prev = orgData(idx).BSgene(i).n(aa)-orgData(idx).p(aa)*orgData(idx).l_avg;
            orgData(idx).BSgene(i).CB(aa) = CB;
            orgData(idx).BSgene(i).CB_prev(aa) = CB_prev;
            
            % choose the minimum for orgData CB and CB_prev
            if CB < orgData(idx).CB(aa)
                orgData(idx).CB(aa) = CB;
            end
            if CB_prev < orgData(idx).CB_prev(aa)
                orgData(idx).CB_prev(aa) = CB_prev;
            end
        end
    end
    % for the rest data, mark as NaN
    for j = 1:20
        if orgData(idx).CB(j) == Inf
            orgData(idx).CB(j) = NaN;
        end
        if orgData(idx).CB_prev(j) == Inf
            orgData(idx).CB_prev(j) = NaN;
        end
    end
end

% organism
% bodysite
% CB -- n-m
% LICB (length independent CB) -- n-pl
% amino acid synthesized
% n -- number of cognate amino acid in the biosyn enzyme
% m -- average number of cognate amino acid in the protein
% p -- percentage of the protein that is the cognate average
% l -- average protein length