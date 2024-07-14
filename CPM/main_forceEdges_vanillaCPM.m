%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  N is the number of subjects
%  T is the number of tasks
%  x: 268*268*N*T
%  y: N*1
%  d: number of elements in lower triangular part of 268*268
%  group: 1*N of type subject
%  subject: a class of 9 task-based connectome and single label
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





function [results] = main_forceEdges_vanillaCPM(x,y,k,edges_pos,edges_neg,dset,mask)
    
    dataset = dset; 
    g = buildGroup(x,dataset,'none'); % mask=false, Bins
    options = [];
    options.thresh=0.05;
    options.seed = randi([1 10000]);
    options.k = k;
    options.control=-1;
    options.phenotype = phenotype('behav',y);
    options.diagnosis = randi(2,175,1);
    options.pos_edges = edges_pos;
    options.neg_edges = edges_neg;

    if mask == 1
        all_mask = find(edges_pos == 0 & edges_neg == 0);
        g.all_edges(all_mask,:) = [];
        g.num_edge = size(g.all_edges,1);
    end

    m = cpm_forceEdges(g,options);
    m.run();
    results = m.evaluate();
end

function g = buildGroup(x,dataset,mask)
N =size(x,3);
subjects(1,N) = subject(N);
for i=1:N
    subjects(i) = subject(x(:,:,i,:),i,dataset,mask);
end
g = group(subjects);
end