function imagequick(x,y,F,Flims)
% imagequick(x,y,F,Flims);
% Uses imagesc but adjusts axes and scales correctly.  Same calling args
% as imagesc.

if (nargin <= 3)
    Flims = [-40 0] + max(real(F(:)));
end

colormap('jet');
imagesc(x,y,F,Flims);

set(gca,'YDir','normal');
colorbar('EastOutside');
axis normal

end