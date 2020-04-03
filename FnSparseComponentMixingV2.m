% Sparse Componnet Mixing Function: FnSparseComponentMixing
%% Inputs:
% m: The number of sensors a.k.a the ambient dimension
% n: The number of sources
% k: The number of acitive source in each time point (nonzero elements
% in each col of S) or the dimension of each subspace
% N: The cardinal of each subspace in kSCA mode
% Nk: the cardinal of  each subspace in MSCA mode
% MixingMode: Sparse Component Mixing Mode:'kSCA'=> Uniform Sparse Component
% Mixing,'MSCA'=> Multiple Sparse Component Mixing
%% Outputs:
% X: The mixtures / mixed data/observed data/data points : m*T
% S: The source matrix/The sparse componenets/The weights in linaer
% combination of the columns of A a basis of the subspaces:n*T
% A: The mixing matrix/ the basis of subspaces/ The dictionary: m*n
% OrthA: The orthonormalized A: each column of OrthA is obtaind of
% orth(A(:,P(i,:)) which P(i,:) is the inices of selected column of A:
% Labels: ground-Truth for subspace clustering:1*T,Label(i) is a number
% between 1 and the total number of the subspaces
% SubspaceInds: each row of P shows the indices of the columns of A that
% has generated each subspace.P is c*k in kSCA mode which c is the total
% number of subspace and is a cell in MSCA mode.
% SubspaceNum: The total number of the subspaces.
% *************************************************************************
% Ehsan Eqlimi, 26 Bahman 1394/ Jan 2015
% Edit 1:Add a parameter to control the variance of gaussian noise over the
% inactive(nonzeros) sources: 9 Xordad 1394 and Add PermkSCA and kSCANoisy
% *************************************************************************


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

function [X,S,OrthA,A,Labels,SubspaceInds,SubspaceNum]=FnSparseComponentMixing(m,n,k,T,N,Nk,Var,IterNum,RankTh,AMode,MixingMode)
if AMode
BestKa=0;
Iter=0;
while BestKa<RankTh && Iter<IterNum
    A=FnColNormalizer(randn(m,n));
    [R,Ka]=k_rank_EE(A);
    if Ka>BestKa
        BestKa=Ka;
        BestA=A;
    end
    Iter=Iter+1;
end
A=BestA;
else
A=FnColNormalizer(randn(m,n));
end
% A=randn(m,n);
% A=FnColNormalizer(A);
P=nchoosek(1:n,k);
c=length(P);
OrthA=[];
% S=[];
X=[];
Labels=[];
S=[];
switch MixingMode
    case 'kSCA'
        SubspaceInds=P;
        for i=1:c
            %             a=orth(A(:,P(i,:)));
            a=(A(:,P(i,:)));
            OrthA=[OrthA a];
            s=randn(k,N(i));
            sWithZeros=zeros(n,N(i));
            sWithZeros(P(i,:),:)=s;
            IdxActive = zeros(n,1);
            IdxActive(P(i,:))=1;
            IdxActive=logical(IdxActive);
            sWithZeros(~IdxActive,:)=Var*randn(n-k,1);
            S=[S sWithZeros];
            x=a*s;
            X=[X x];
            %             S=[S s];
            Labels=[Labels i*ones(1,N(i))];
        end
        
    case 'PermkSCA'
        SubspaceInds=P;
        S=zeros(n,T);
        for i=1:T
            IdxActive = zeros(n,1);
            %         IdxRandPerm = randperm(k);
            Labels(i) = randperm(size(P,1),1);
            IdxRandPerm=P(Labels(i),:);
            IdxActive(IdxRandPerm) = 1;
            IdxActive = logical(IdxActive);
            S(IdxActive,i) = randn(k,1);
            
            S(~IdxActive,i) = Var*randn(n-k,1);
        end
        X=A*S;
    case 'kSCANoisy'
      SubspaceInds=P;
        for i=1:c
            %             a=orth(A(:,P(i,:)));
% %             a=(A(:,P(i,:)));
% %             OrthA=[OrthA a];
            s=randn(k,N(i));
            sWithZeros=zeros(n,N(i));
            sWithZeros(P(i,:),:)=s;
            IdxActive = zeros(n,1);
            IdxActive(P(i,:))=1;
            IdxActive=logical(IdxActive);
            sWithZeros(~IdxActive,:)=Var*randn(n-k,1);
            S=[S sWithZeros];
%             x=a*s;
%             X=[X x];
            %             S=[S s];
            Labels=[Labels i*ones(1,N(i))];
        end
        X=A*S;
    case 'MSCA'
        Labels=0;
        for j=1:k
            Added(j)=max(Labels);
            kj=j;
            Pj=nchoosek(1:n,kj);
            SubspaceInds{j}=Pj;
            cj=length(Pj);
            for i=1:cj
                a=orth(A(:,Pj(i,:)));
                OrthA=[OrthA a];
                s=randn(kj,Nk(i,j));
                sWithZeros=zeros(n,Nk(i,j));
                sWithZeros(Pj(i,:),:)=s;
                
                S=[S sWithZeros];
                x=a*s;
                X=[X x];
                %                 S=[S s];
                Labels=[Labels (i+Added(j))*ones(1,Nk(i,j))];
                
            end
            
        end
        Labels=Labels(2:end);
end
SubspaceNum=max(Labels);