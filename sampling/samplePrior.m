function VARSlist = samplePrior(PRIOR,n)
names = fieldnames(PRIOR);
nparams = length(names);
VARSlist = cell(n,1);

%Initialize
for i=1:n
    for k=1:nparams
        VARSlist{i} = setfield(VARSlist{i},names{k},nan);
    end
end

%Populate
for k = 1:nparams
    %See if value of field is numeric and treat as a constant,
    val = getfield(PRIOR,names{k});
    if isa(val,'numeric')
        for i=1:n
            VARSlist{i} = setfield(VARSlist{i},names{k},val);
        end
        
    %Otherwise assume its a probability distribution and sample
    else
        for i=1:n
            VARSlist{i} = setfield(VARSlist{i},names{k},random(getfield(PRIOR,names{k})));
        end
    end
end
