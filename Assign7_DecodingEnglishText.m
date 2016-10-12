%% code snippets

% snippet for converting text into a double array
vtext = double(text) - 96;
index = find(vtext==-64);
if (~isempty(index))
    vtext(index)=27;
end

% convert the message in numbers to text
ctext = cvtext;
index = find(ctext==27);
if (~isempty(index))
    ctext(index)=-64;
end
ctext = char(ctext+96);


% load the English transition matrices and cleaning them
load English_trans;

% A = bigram transition A(i,j) = Prob(x_t = i| x_{t-1} = j)
% Clean A -- remove zeros from the transition matrix and renormalize
tol = 1e-6;
A = A + tol;
for k = 1:size(A,2)
    A(:,k) = A(:,k)/sum(A(:,k));
end


% S = trigram transition A(i,j,k) = Prob(x_t = i| x_{t-1} = j, x_{t-2} = k)
% Clean S -- remove zeros from the transition matrix and renormalize
S = S + tol;
for k = 1:27
    for l = 1:27
        sumprob = 0;
        for m = 1:27
            sumprob = sumprob + S(m,l,k);
        end
        if (sumprob>0)
            for m = 1:27
                S(m,l,k) = S(m,l,k)/sumprob;
            end
        end
    end
end
    

% compute the log-likelihood using only bigrams A
% given a permutation xterm and coded text cvtext
llxperm = 
for k = 2:length(cvtext)
    llxperm = llxperm + log(A(xperm(cvtext(k)),xperm(cvtext(k-1))));
end


% compute the log-likelihood sing only trigrams S
% given a permutation xterm and coded text cvtext
llxperm = 0; 
for k = 3:length(cvtext)
    llxperm = llxperm + log(S(xperm(cvtext(k)),xperm(cvtext(k-1)),xperm(cvtext(k-2))));
end

%% coded text that you need to decode

ctext = 'sopxrqprphastzopfinbpbrupscprmasnprcbpoz pfnifjqpx a pqoasjsctpozsao  c'
cvtext = [19
    15
    16
    24
    18
    17
    16
    18
    16
    8
    1
    19
    20
    26
    15
    16
     6
     9
    14
     2
    16
     2
    18
    21
    16
    19
     3
    16
    18
    13
     1
    19
    14
    16
    18
     3
     2
    16
    15
    26
    27
    16
     6
    14
     9
     6
    10
    17
    16
    24
    27
     1
    27
    16
    17
    15
     1
    19
    10
    19
     3
    20
    16
    15
    26
    19
     1
    15
    27
    27
     3]';

