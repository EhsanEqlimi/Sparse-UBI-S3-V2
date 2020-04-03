
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

function [X,SWithZeros,OrthA,A,Labels]=FnSparseComponentMixing(m,n,k,N,Nk,Mode)
A=FnNormalizer(randn(m,n));
P=nchoosek(1:n,k);
c=length(P);
OrthA=[];
% S=[];
X=[];
Labels=[];
SWithZeros=[];
switch Mode
    case 'kSCA'
        for i=1:c
            a=orth(A(:,P(i,:)));
            OrthA=[OrthA a];
            s=randn(k,N(i));
            sWithZeros=zeros(n,N(i));
            sWithZeros(P(i,:),:)=s;
            SWithZeros=[SWithZeros sWithZeros];
            x=a*s;
            X=[X x];
%             S=[S s];
            Labels=[Labels i*ones(1,N(i))];
        end
    case 'MSCA'
        Labels=0;
        for j=1:k
            Added(j)=max(Labels);
            kj=j;
            Pj=nchoosek(1:n,kj);
            cj=length(Pj);
            for i=1:cj
                a=orth(A(:,Pj(i,:)));
                OrthA=[OrthA a];
                s=randn(kj,Nk(i,j));
                sWithZeros=zeros(n,Nk(i,j));
                sWithZeros(Pj(i,:),:)=s;
                SWithZeros=[SWithZeros sWithZeros];
                x=a*s;
                X=[X x];
                %                 S=[S s];
                Labels=[Labels (i+Added(j))*ones(1,Nk(i,j))];
                
            end
            
        end
        Labels=Labels(2:end);
end