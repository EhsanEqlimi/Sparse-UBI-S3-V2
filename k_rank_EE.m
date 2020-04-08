%Written by Dr.Bahador Makkiabadi, @TUMS,tehran, Iran
% [1] E. Eqlimi and B. Makkiabadi, “Multiple sparse component analysis
% based on subspace selective search algorithm,” in 2015 23rd Iranian
% Conference on Electrical Engineering. IEEE, 2015, pp. 550–554
% 
% [2] E. Eqlimi and B. Makkiabadi, “An efficient K-SCA based unerdetermined
% channel identification algorithm for online applications,” in 2015
% 23rd European Signal Processing Conference (EUSIPCO). IEEE, 2015,
% pp. 2661–2665.

%Please cite the above papers in case of any usage or benchmarking.
% 
% (C) Bahador Makkiabdi, Ehsan Eqlimi, Jan 2015, @Biomedical Engineering Group, Tehran
% Universty of Medical Scineces (TUMS), Tehran, Iran
% https://github.com/EhsanEqlimi/Sparse-UBI-S3

% Ehsun.Eghlimi@gmail.com, Eghlimi@razi.tumss.ac.ir
%% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
% OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
% AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
% IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
% OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [min_krank,mineig,min_sel,minvec,all_eig,all_sel,all_vec,worst_case,krank_old,mineig_old,min_sel_old,minvec_old]=k_rank_EE(A,thr,maxlevel)
n=size(A,2);
m=size(A,1);

if(~exist('maxlevel'))
   maxlevel=min(n,m); 
end

%min_krank=maxlevel;
min_krank=min(m,n);
worst_case=[];

if(~exist('thr'))
thr=1e-3;
end

mineig=1e20;
minvec=[];
min_sel=[];
krank=1;
st=1;
for i=1:min(n,min(m,maxlevel))%min(n,m):-1:2

           if(st)
               ii=min(n,min(m,maxlevel));
%          out=combntns(1:n,ii); 
           out=nchoosek(1:n,ii);
         j=1;
                svx=mnorm(svd(A(:,out(j,:)))) ;

       all_eig(:,j,ii)=mineig;    
       
       all_vec(:,j,ii)=svx; 
       all_sel(:,j,ii)=out(j,:);    
       st=0;
           end
    
        
%if(m<=n)
% out=combntns(1:n,i);
out=nchoosek(1:n,i);
%end
    for j=1:size(out,1)
       svx=mnorm(svd(A(:,out(j,:)))) ;
       if(svx(end)<mineig)
           
          mineig_old=mineig;
          minvec_old=minvec;
          min_sel_old=min_sel;
          krank_old=krank;
           
          mineig=svx(end);
          minvec=svx;
          min_sel=out(j,:);
          krank=i;
          
       end
       
       all_sel(1:size(out,2),j,i)=out(j,:);    
       all_vec(1:length(svx),j,i)=svx;
       all_eig(:,j,i)=svx(end);    
       
if(svx(end)<=thr),
    min_krank=min(i-1,min_krank);
    zx=zeros(1,maxlevel);
    zx(1:length(out(j,:)))=out(j,:);
    worst_case=[worst_case;zx];
    if(thr>0)
    break;
    end
    
end

           end
    
if(mineig<=thr),
   % worst_case=out;
   if(thr>0)
    break;
   end
end

end
%krank,mineig,min_sel,minvec
if 0&(~isempty(worst_case))
           mineig=mineig_old;
          minvec=minvec_old;
          min_sel=min_sel_old;
          krank=i-1;%krank_old;
end

if (~isempty(worst_case))
%krank=i-1;
end
    
    
