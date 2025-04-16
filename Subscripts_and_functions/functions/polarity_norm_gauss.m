%function to normalise polarity of a time series of gauss coefficient

function gh_full_n = polarity_norm_gauss(gh_full)

Ntimesteps = size(gh_full,2);
gh_full_n=zeros(size(gh_full)); gh_full_n =[];


for i=1:Ntimesteps
    if(gh_full(1,i)) > 0
        gh_full_n(:,i)=gh_full(:,i).*-1;
    else
        gh_full_n(:,i)=gh_full(:,i);
    end
end