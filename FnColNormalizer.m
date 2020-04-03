% Normalize every cloumn of X to [0,1] range
%==========================================================================
% Syntax: X=FnNormalizer(X)
%==========================================================================
% Input variable: 
% 1)X : Input Matrix or column vector.The columns are observation and the
% rows are 
%==========================================================================
% Output variable:
% 1)X: The matrix after normalization
%==========================================================================
% Ehsan Eqlimi (Ehsun.Eghlimi@gmail.com)
% Date: 20-October-2014
%==========================================================================

% [1] E. Eqlimi and B. Makkiabadi, �Multiple sparse component analysis
% based on subspace selective search algorithm,� in 2015 23rd Iranian
% Conference on Electrical Engineering. IEEE, 2015, pp. 550�554
% 
% [2] E. Eqlimi and B. Makkiabadi, �An efficient K-SCA based unerdetermined
% channel identification algorithm for online applications,� in 2015
% 23rd European Signal Processing Conference (EUSIPCO). IEEE, 2015,
% pp. 2661�2665.

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

function X=FnColNormalizer(X)

for i=1:size(X,2)
    X(:,i)=X(:,i)/norm(X(:,i));
end