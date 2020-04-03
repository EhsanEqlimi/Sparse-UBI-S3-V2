
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

function [MixingIdentificatioError, MixingVectorerror,NMSE,NMSSum,AhatNew,Norm2Error] = FnMixingIdentificationError(A,Ahat)
%% Error estimation
[AhatNew,M]=nearest2(Ahat,(A));
BAS=0;
n=min([size(A,2) ,size(AhatNew,2)]);
for e=1:n
    MixingVectorerror(e)=acos((AhatNew(:,e)'* A(:,e))/(norm(AhatNew(:,e))*norm(A(:,e))));
    BAS=MixingVectorerror(e)+BAS;
end
MixingIdentificatioError=abs((180/pi)*(BAS))
% AhatNew
% AhatNew

%% NMSE
[NMSE,NMSSum]=FnNMSECalc(A,AhatNew);
NMSSum
%% norm 2 Error
Norm2Error=norm(A-AhatNew)
