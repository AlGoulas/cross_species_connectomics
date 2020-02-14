


function w_r=Rescale01(w)

    maxx    = max(w);
    minxx   = min(w);
    w_r  = (w - repmat(minxx,size(w,1),1)) ./ repmat(maxx-minxx,size(w,1),1);


return