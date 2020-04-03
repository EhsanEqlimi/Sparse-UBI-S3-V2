
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
% (C) Ehsan Eqlimi, Jan 2015, @Biomedical Engineering Group, Tehran
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

function [Ahat]=FnCS2(C,m,Ahat,Th)
F=0;
if size(C,2)> m-1
    AllPerms=nchoosek(1:size(C,2),m);
    for i=1:size(AllPerms,1)
        Selected = C(:,AllPerms(i,:));
        [V,D]=eig(Selected*Selected');
        v=nchoosek(1:m,m-1);
        vv=(V(:,1));
        d=diag(D);
        d=d(2)/d(end);
        d3=max(max(abs(vv'*Selected)));
        for j=1:size(v,1)
            [V1,temp]=eig(Selected(:,v(j,:))*Selected(:,v(j,:))');
            temp=abs(diag(temp));
            d(j+1)=temp(2)/temp(end);
            vv(:,j+1)=V1(:,1);
            d3(j+1)=max(max(abs(V1(:,1)'*Selected)));
        end
        D=diag(D);
        D=abs(D);
        if abs(D(1)/D(end))<Th  && abs(D(2)/(D(end)+eps))>Th*1e8 %abs(D(2)/(D(1)+D(2)))>.999
            if max(d3)<Th %min(abs(d))>Th*1e8 && dd(end-1)/dd(end)<Th*1e4
                F=F+1;
                Ahat(:,F)=V(:,1);
                
            end
        end
    end
end