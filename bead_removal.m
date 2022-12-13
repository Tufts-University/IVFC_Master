function [M] = bead_removal(i, M, beadfile)
    load(beadfile,'peak_values')
    locs = peak_values(peak_values(:,6)==i,7);
    width = round(peak_values(peak_values(:,6)==i,17).*1.5);
    means = mean(M,1);
    for l = 1:length(locs)
        M(locs(l)-width(l):locs(l)+width(l),:) = repmat(means,width(l)*2+1,1);
    end
end