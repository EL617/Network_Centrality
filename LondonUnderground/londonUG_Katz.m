%LONDON_UG.M
%
% LOAD LONDON TRANSPORT DATA
% FROM
% Nonlocal PageRank
% Tudisco, Cipolla and Durastante, 
% to appear in  
% ESAIM: Mathematical Modelling and Numerical Analysis
% https://arxiv.org/abs/2001.10421
% Github link for data is at
% https://github.com/Cirdans-Home/NonLocalPageRank
%
% This code is written to be as simple and understandable
% as possible. Not made for efficiency or elegance.
%
% DJH Nov 2020
%


% Load the data
% Note: initial line of descriptive characters 
% has been removed from these files
%
% Further comments:
%There are 369 (nodes) stations in the edge data,
%/london_tube_edges.txt
%
% labelled from 0 to 368.
% But only 271 stations have corresponding 
% passenger data---these are for the undergound 
% train lines.
%
% These 271 stations appear in the passenger usage list at
%/london_tube_usage.txt
%

Edges = load('cropped_edges.txt');
%%% First row of Edges is [4     1    77     1]
% 4 is id of the transport line 
% (we need this to be between 1 and 11)
% edge is from node 1 to node 77
% weight is 1 (as for all edges here)

% Create adjacency matrix 
% We know that there are 369 nodes
% ids go from 0 to 368
% but MATLAB starts indices at 1 not 0
Nbig = 369;
Abig = zeros(Nbig,Nbig);
for count = 1:length(Edges) % go over all rows
    if Edges(count,1) < 12      % underground links only 
                                           % change if we want all links
        ival = Edges(count,2);
        jval = Edges(count,3);
       Abig(ival+1,jval+1) = 1;    % insert the edge
    end
end
Abig = sign(Abig+Abig');       % make edges symmetric
% Abig is the adjacency matrix for the network.
% Note that Abig has an empty row and column for
% any station that is not part of the underground network.
% (That is, for any link that had transport line 12 or 13.)
% There are 98 of these rows: sum(sum(Abig)==0) = 98
% This is OK, as the Katz centrality for these nodes does not
% affect the others.

% Now look at usage data
% We only have usage data for the underground stations,
% these are 271 nodes
Usage = importdata('cropped_usage.txt');
Unum = Usage.data;
%%% First row of Unum is 
% 6.0405    6.2746    6.2350    6.0590    5.7490    5.5849...
% 5.4624    5.4283    5.4595    5.8168
% These are passenger numbers, in millions, for the years
% 2017 2016 2015 2014 2013 2012 2011 2010 2009 2008
%
% Usage.textdata is a 271×2 cell array
% with first "row"
%   {'4'  }    {'actontown'                 }
% We need to pick out the id 4 as an integer
for k = 1:length(Usage.textdata)
     stid(k) = str2num(Usage.textdata{k});
end
% First entry of stid is 4 (the 5th station).
% This is Acton Town.
% It is the station id for first row of Unum

% For later use, we can create a list of the rows and 
% columns in Abig that are completely zero:
% i.e., correspond to non-underground stations.
% We will do this by finding all entries missing from stid
idmiss = [];
for k = 1:Nbig
    if ~ismember(k,stid+1)
        idmiss = [idmiss,k]; % add k to this list of ids
    end
end

%sanity check
% length(idmiss) gives 98
%max(max(Abig(idmiss,idmiss))) gives 0

% max alpha = 0.2644;
alpha = 0.25;
% Set the comparison range;
range = 100;

% Find the Katz Component for each station
Katz_component = (eye(Nbig,Nbig) - alpha*Abig)\(ones(Nbig,1));
% Here using the average number of passengers per year
avg_num = mean(Unum,2);
% This gives us a 271 x 1 matrix of average number of passengers in the
% order of stid

% Set range of our measure (First k max numbers)
k_max = 271;
% Find the station id of the maximum Katz Component from the matrix
[max_katz,kz_idxs] = maxk(Katz_component,k_max);
% Use the index in the matrix to find the corresponding index from stid
% which will be the same index for the station id in avg_num
max_st_kz = [];
for i = 1:k_max
    idx = find(stid == (kz_idxs(i) - 1));
    max_st_kz(i,1) = idx;
end
% This would be the station id with the k_max highest Katz centrality
% measure

% Find the average number of passengers from the stations
p_num_kz = [];
for i = 1:k_max
    p_num_kz(1,i) = avg_num(max_st_kz(i),1);
end
% This is the number of passengers from k_max highest Katz centrality
% measure

% Normalise the data
katz_data = normalize(max_katz);
p_num_data_kz = normalize(p_num_kz');

% Find correlation
% Pearson Correlation
r_ps_mat_kz = corrcoef(p_num_data_kz(1:range,:),katz_data(1:range,:));
r_ps_kz = r_ps_mat_kz(1,2);
% Spearman Correlation
r_sp_kz = corr(p_num_data_kz(1:range,:),katz_data(1:range,:),'type','Spearman');
% Kendall Correlation
r_kn_kz = corr(p_num_data_kz(1:range,:),katz_data(1:range,:),'type','Kendall');

% Here looking for the overlap of stations occur in the range of k_max
% Represent as a ratio
% First find the station id of the maximum amount of passengers
[max_pas_num, p_idxs] = maxk(avg_num,k_max);
max_st_pnum = [];
for i = 1:k_max
   max_st_pnum(i,1) = p_idxs(i);
end
% Calculate the overlap ratio
overlap_kz = intersect(max_st_pnum(1:range,:), max_st_kz(1:range,:));
ol_ratio_kz = length(overlap_kz)/range;

% figure(1)
% clf
% spy(Abig)   % shows nonzero pattern of matrix

% Cumulative Graph
cum_pnum = cumsum(max_pas_num);
cum_katz = cumsum(p_num_kz);



% Try Eigenvector Centrality Measure
sparse_A_big = sparse(Abig);
[e_vecs,e_vals] = eigs(sparse_A_big);
vals = transpose(diag(e_vals));
[max_val,idx] = max(real(vals));
eigen_component = abs(real(e_vecs(:,idx)));
% Find the station id of the maximum Eigenvector Component from the matrix
[max_eigen,eg_idxs] = maxk(eigen_component,k_max);
% Use the index in the matrix to find the corresponding index from stid
% which will be the same index for the station id in avg_num
max_st_eg = [];
for i = 1:k_max
    idx = find(stid == (eg_idxs(i) - 1));
    max_st_eg(i,1) = idx;
end
% This would be the station id with the k_max highest Eigenvector centrality
% measure

% Find the average number of passengers from the stations
p_num_eg = [];
for i = 1:k_max
    p_num_eg(1,i) = avg_num(max_st_eg(i),1);
end
% This is the number of passengers from k_max highest Eigenvector centrality
% measure
eigen_data = normalize(max_eigen);
p_num_data_eg = normalize(p_num_eg');

% Find correlation
% Pearson Correlation
r_ps_mat_eg = corrcoef(p_num_data_eg(1:range,:),eigen_data(1:range,:));
r_ps_eg = r_ps_mat_eg(1,2);
% Spearman Correlation
r_sp_eg = corr(p_num_data_eg(1:range,:),eigen_data(1:range,:),'type','Spearman');
% Kendall Correlation
r_kn_eg = corr(p_num_data_eg(1:range,:),eigen_data(1:range,:),'type','Kendall');


% Calculate the overlap ratio
overlap_eg = intersect(max_st_pnum(1:range,:), max_st_eg(1:range,:));
ol_ratio_eg = length(overlap_eg)/range;

% Cumulative Graph
cum_eigen = cumsum(p_num_eg);



% Try the Degree centrality of each station
degree_component = sum(Abig,2);
% Find the station id of the maximum Degree Component from the matrix
[max_degree,dg_idxs] = maxk(degree_component,k_max);
% Use the index in the matrix to find the corresponding index from stid
% which will be the same index for the station id in avg_num
max_st_dg = [];
for i = 1:k_max
    idx = find(stid == (dg_idxs(i) - 1));
    max_st_dg(i,1) = idx;
end
% This would be the station id with the k_max highest Degree centrality
% measure
% Find the average number of passengers from the stations
p_num_dg = [];
for i = 1:k_max
    p_num_dg(1,i) = avg_num(max_st_dg(i),1);
end
% This is the number of passengers from k_max highest Eigenvector centrality
% measure
% Normalise the data
degree_data = normalize(max_degree);
p_num_data_dg = normalize(p_num_dg');
% Cumulative Graph
cum_degree = cumsum(p_num_dg);

% Calculate the overlap ratio
overlap_dg = intersect(max_st_pnum(1:range,:), max_st_dg(1:range,:));
ol_ratio_dg = length(overlap_dg)/range;

% Find correlation
% Pearson Correlation
r_ps_mat_dg = corrcoef(p_num_data_dg(1:range,:),degree_data(1:range,:));
r_ps_dg = r_ps_mat_dg(1,2);
% Spearman Correlation
r_sp_dg = corr(p_num_data_dg(1:range,:),degree_data(1:range,:),'type','Spearman');
% Kendall Correlation
r_kn_dg = corr(p_num_data_dg(1:range,:),degree_data(1:range,:),'type','Kendall');

figure(2)
clf
Katz_cent = Katz_component(stid,:);
plot(stid,Katz_cent,'*')
xlabel('Station Index') 
ylabel('Katz Centrality') 
title('Katz Centrality')

figure(3)
clf
eigen_cent = eigen_component(stid,:);
plot(stid,eigen_cent,'*')
xlabel('Station Index') 
ylabel('Eigenvector Centrality') 
title('Eigenvector Centrality')

figure(4)
clf
degree_cent = degree_component(stid,:);
plot(stid,degree_cent,'*')
xlabel('Station Index') 
ylabel('Degree Centrality') 
title('Degree Centrality')




figure(6)
clf
scatter(p_num_data_kz(1:range,:),katz_data(1:range,:))


figure(7)
clf
plot(cum_pnum,'k')
hold on
plot(cum_eigen,'b')
hold on
plot(cum_katz,'g')
hold on
plot(cum_degree,'r')
legend({'# of Passenger','Eigenvector','Katz','Degree'},'Location','east');
xlabel('Stations') 
ylabel('Cumulative Passenger Number in Millions') 
title('Cumulative Comparison')


% figure(6)
% clf
% surf(cum_eigen,cum_katz,cum_pnum_eg')
% title('Passengers vs Centrality measures')
% xlabel('index')
% ylabel('#')
% set(gca,'FontWeight','Bold')
% set(gca,'FontSize',18)

% ax = gca;
% Requires R2020a or later
% exportgraphics(ax,'iprplot.png','Resolution',300) 

