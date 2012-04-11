function ltid_chk_vwg(v,w,g,o)
% function ltid_chk_vwg(v,w,g,o)
%
% visual check of quality of matching v,w data by g 
% with o!=0 displays real parts only

if nargin<3, error('3 inputs required'); end
if nargin<4, o=0; end
for i=1:size(v,2),
    if o==0,
        close(gcf);
        subplot(2,1,1);plot(w,real(v(:,i)),'.',w,real(g(:,i))); grid
        subplot(2,1,2);plot(w,imag(v(:,i)),'.',w,imag(g(:,i))); grid
    else
        close(gcf);plot(w,real(v(:,i)),'.',w,real(g(:,i))); grid
    end
    x=input('Press ENTER to continue ... ');
    if ~isempty(x), break; end
end
