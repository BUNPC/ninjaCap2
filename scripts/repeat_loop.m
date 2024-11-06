function [ans,hHex_outline_adjust] = repeat_loop(vHex2,eHex,hHex,hHex_outline_adjust,eHexLen,vHex_outline_idx)
     outline_lengths_init = [];
     outline_lengths = [];
     final_idx = zeros(size(eHex,1),1);
     for u = 1:length(vHex_outline_idx)-1
        e =  [vHex_outline_idx(u) vHex_outline_idx(u+1)];
        idx = ismember(eHex,e,'row');
        final_idx = final_idx | idx;
     end
     outline_lengths_init = hHex(final_idx);
    outline_lengths = eHexLen(final_idx);
     length_err = outline_lengths-outline_lengths_init;
     if sum(length_err) > 2 || any(length_err>0.2)
         ans = true;
         hHex_outline_adjust(final_idx) = hHex_outline_adjust(final_idx)-length_err;
     else
         ans = false;
     end
end