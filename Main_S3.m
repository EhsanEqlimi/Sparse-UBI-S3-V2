%% Main cos for running S3 method 
% this cose generates simulated mixtures (X=A.S), and then applies
% S3 method on it, in order to identify the mixing matrix (A). It also calculate 
% the idenrtification error.

%This code implemntes our underderminded blind identification method 
%proposed in the follwong papers:

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

clc
clear
close all
%% Preliminaries
m=4; % The number of sensors a.k.a the ambient dimension
n=5; % The number of sources a.k.a
k=m-1; % The number of acitive source in each time point (nonzero elements
%in each col of S) or the dimension of each subspace
r=m-k; % The dimension of the normal space (normal vector/orthogonal complement) of each subspace
c=nchoosek(n,k); % The number of all possible r-dim subspaces
T= 2000; % The Number of data poins a.k.a the number of samples (time points)
Th=1e-2; % Remove close to origin
ReNum=300;% The number of recieved data
Th1=1e-15; % for RANSAC-Like Subspace Identifying
Th2=1e-5; % for Clustering on Detected Subspaces based on Norm not EVD
Th3=1e-14;
% FuzzyArt Parameters
px=.6;
alphax=.09;
hypin=1;
Sigma=0; % Parameter controls the standard deviation of normal noise over the
% zero sources (non-strict saprsity)
AMode=1; % k-SCA Condition for A are satisfied
IterNum=500;% The number of iteration to generate a good mixing matrix
RankTh=0.1; % to generate a good mixing matrix
% Sigma1=0.02;
MixingMode='kSCA'%'MSCA'%'kSCA';%'kSCA';% Uniform Sparse Component Mixing  %'MSCA'=>Multiple Sparse Component Mixing
N=zeros(1,c);
N(1,:)=ceil(T/c);
Nk=zeros(c,k);
for j=1:k
    for i=1:nchoosek(n,j)
        Nk(i,j)=ceil(ceil(T/k)/nchoosek(n,j));% N showes the number of each subspace
    end
end

%% Mixing
%% SCA Mixing (Simulted Sources)
%%'kSCA';%'kSCA';% Uniform Sparse Component Mixing  %'MSCA'=>Multiple Sparse Component Mixing, 'PermKSCA'
MixingMode='PermkSCA'% {'kSCANoisy'%'PermkSCA'%kSCA';%'MSCA'};
Orth=0; %if Orth=1 A is orth
N=zeros(1,c);
N(1,:)=ceil(T/c); % N is the number of the subspaces in k-SCA mode
Nk=zeros(c,k);
for j=1:k
    for i=1:nchoosek(n,j)
        Nk(i,j)=ceil(ceil(T/k)/nchoosek(n,j));% N showes the number of each subspace in MSCA mode
    end
end

%% Sparse Component Mixing
[X,S,OrthA,A,Labels,SubspaceInds,SubspaceNum]=FnSparseComponentMixingV2(m,n,k,T,N,Nk,Sigma,IterNum,RankTh,AMode,MixingMode);
% [X,SWithZeros,OrthA,A,Labels]=FnSparseComponentMixing(m,n,k,N,Nk,MixingMode);
%% Remove the cose to origin data points
X=FnRemoveColCloseOrigin( X,Th );
%% Normalization to norm
X=FnNormalizer(X);
NormVecSubspaces=[];
NMax=60;
CountinEachWhile=0;
Level=0;
CentroidNormVecSubpaces=[];
Iter=0;
Ahat=[];
C_Ahat=[];
AllAhat=[];
N=0;
while size(C_Ahat)<n %N<NMax 
    
    [NormVecSubspaces,CountinEachWhile,NewNormVecSubspace]=FnS3(X,Th1,Level,NormVecSubspaces,CountinEachWhile); %Alg1:S3
    
    if CountinEachWhile~=0;
        CountinEachWhile=0;
        Iter=Iter+1;
        [CentroidNormVecSubpaces,Winner]=FnBBC(NewNormVecSubspace,CentroidNormVecSubpaces,Iter,Th2); %Alg2:BBC
        [Ahat]=FnCS2(CentroidNormVecSubpaces,m,Ahat,Th3);%Alg3:CS2
        if size(Ahat,2)~=0
            N=N+1;
            for j=1:size(Ahat,2)
                [C_Ahat,Winner]=FnBBC(Ahat(:,j),C_Ahat,N,Th1);%Alg2:BBC
                
                ChannelNum=size(C_Ahat,2)
               
                AllAhat=[AllAhat Ahat];
            end
        end
    end
end
[MixingIdentificatioError, MixingVectorerror,NMSE,NMSSum,AhatNew,Norm2Error]=FnMixingIdentificationError(A,C_Ahat);

A
AhatNew
