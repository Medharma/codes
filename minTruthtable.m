%function [Bins,inps,Nums,ott] = minTruthtable(tt, flags)
% Function to minimize a truth table using the Quine-McCluskey algorithm.
% * tt = truth table, as character (filled with '0's, '-'s or '1's). E.g. '00011-10'.
%        tt must not be more than 2^15 long, due to an int16 optimization.
% * flags = character array with the following control flags:
%   - 'e': do an exact optimization, which will guarantee the optimum (with
%          reservation for potential bugs in the code).
%   - 'Q': do a quick heuristic, as defined by the Quine-McCluskey, which is
%          much faster for big problems, but does not guarantee the optimum.
%   - 'q': do another quick heuristic.
%   - '<number>': A semi-quick heuristic. '0': fastest, higher: better result
%          (run <number> levels of binary search, and below this, the 'q' heuristic).
%   - 'v': Verify the result using a self-check (make sure the reported result
%          meets the provided truth table).
%   Backward compatibility: If flags is boolean => see the 'v' flag.
% * Bins = character matrix. E.g. ['1-0'; '011'] (meaning => x2*/x0 + /x2*x1*x0).
% * inps = number of logical gates input for a minimal nand-nand-synthesis of tt. E.g. 7
% * Nums = cell array with arrays indicating which terms contains which '1's in the tt.
%          (0-indexed, so '1001' => [0 3]. E.g. {[4 6], [3]}.
% * ott = output truth table, as character, filled with '0's and '1'. E.g. tt='00011010'
% 
% If none of the optimization flags are provided, a suitable will be used,
% considering problem size.
% 
% Examples:
% # >> minTruthtable('0-11--10110010--')
% # | Kmap (index to left, -:unused don't care, =:used don't care):
% # |  /  0  1  3  2 \  / . - 1 1 \
% # |  |  4  5  7  6 |  | = - . 1 |
% # |  | 12 13 15 14 |  | 1 . - = |
% # |  \  8  9 11 10 /  \ 1 1 . . /
% # | All terms:
% # | * T( 1): "001-" <-> {2 3}
% # | * T( 2): "100-" <-> {8 9}
% # | * T( 3): "-1-0" <-> {6 4 12 14}
% # | y = ~x(3)*~x(2)*x(1) + x(3)*~x(2)*~x(1) + x(2)*~x(0);
% # | 
% # | Logical complexity: 11 inputs
% # | 
% # | Input tt:  '0-11--10110010--'
% # | Output tt: '0011101011001010'
% # >> [Bins,inps,Nums,tt] = minTruthtable('0-11--10110010--');
% -> Bins = ['001-'; '100-'; '-1-0']
% -> inps = 11
% -> Nums = {[2 3]; [8 9]; [6 4 12 14]}
% -> tt = '0011101011001010'
% 
% Author: Petter Källström, petterk@isy.liu.se
% Date: June, 2012
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY. Please use the built-in verifier to make
% sure the output terms corresponds to the input truth table (the 'v' flag).
%
% Acknowledgement: This file is inspired by Andey Popov's function "minBool"
% (http://www.tu-harburg.de/~rtsap/#Programs)

function [Bins,inps,Nums,ott] = minTruthtable(tt, flags)
  % Step description: http://en.wikipedia.org/wiki/Quine%E2%80%93McCluskey_algorithm
  
  %% Initial setup
  if nargin==1 && ~nargout && strcmp(tt,'BIST'), BIST(); return; end
  if nargin < 1
    % Use a test vector if nothing was given.
    tt = '0000100-1-1110-1'; % m(4,8,10,11,12,15) + d(7,9,14)
    flags = true;
    % Kmap:  0 1 3 2  <- x10
    %   / 0  0 0 0 0
    % x32 1  1 0 - 0  Sol: x2*/x1*/x0 + x3*/x2 + x3*/x0 + x3*x1
    %   | 3  1 0 1 -
    %   \ 2  1 - 1 1
    % This test vector will be used as example in the code.
  end
  if 0 < nargin && nargin < 2, flags = ''; end
  doexact = false;
  doquick = false;
  dolevel = NaN;
  doquine = false;
  if ischar(flags)
    doexact = any(flags == 'e'); % complete optimization
    doquick = any(flags == 'q'); % quick optimization
    doquine = any(flags == 'Q'); % original quine-mcCluskey
    dolevel = str2double(flags('0'<=flags & flags<='9')); % NaN if missing
    verify = any(flags == 'v');
  end
  if ~doexact && ~doquick && ~doquine && isnan(dolevel)
    % no directives => guess, based on the input size
    if length(tt) <= 32
      doexact = true;
    elseif length(tt) <= 1024
      dolevel = 8;
    else
      dolevel = 5;
    end
  end
  if ~ischar(tt)
    if islogical(tt)
      tt = char('0'+tt);
    elseif isfloat(tt)
      tt2 = char('-'+zero(size(tt)));
      tt2(tt==0) = '0';
      tt2(tt==1) = '1';
      tt = tt2;
    else
      error('minTruthtable:tt_type', 'tt should be a character array. tt is of type %s.', class(tt))
    end
  end
  if size(tt,2) == 1, tt = tt'; end
  assert(size(tt,1) == 1, 'minTruthtable:multi_output', 'Error in minTruthtable: tt must be a vector.');
  N = round(log2(size(tt,2)));
  assert(size(tt,2) == 2^N, 'minTruthtable:tt_size_mismatch', 'Error in minTruthtable: tt has height %d, which is not a power of two.', size(tt,1));
  assert(all(tt=='-' | tt=='0' | tt=='1'), 'minTruthtable:tt_char', 'Error in minTruthtable: Unexpected character in tt');
  
  % N: Number of variables
  % *On*: everything that has to do with '1' in the truth table
  % *Dc*: everything that has to do with '-' in the truth table
  % n*:   count variable of *
  % An "atomic" term: Like a '1' in the Kmap
  % A term: Like a "ring" in the Kmap
  Ons = find(tt=='1'); % '1's. Ex: [4 8 10 11 12 15] + 1
  Dcs = find(tt=='-'); % '-'s. Ex: [7 9 14] + 1
  nOn = length(Ons);
  nDc = length(Dcs);
  
  %% Step 1: Generate product terms
  % Variables used in the loop:
  % Each row in those stands for a product term. Initially atomic terms.
  % - Bins*: binary (char) representations of the terms. Ex: x3*/x0 = '1--0'
  % - Covs*: true iff the term does not have to be included (again). Cov = Covered (or don't care at input).
  % - Nums*: vector with all atomic terms included in the term. Ex: '-100' => [4 12]
  % - Dons*: position of all '-' as int. Eg. Bins='-1-0' => Dons=[1 0 1 0]*[8;4;2;1] = 10
  %          Dons=-1 => This is a double and should be ignored.
  % - Ones*: position of all '1' as int. Eg. Bins='-1-0' => Dons=[0 1 0 0]*[8;4;2;1] = 4
  % - nVars*: Number of rows.
  % - *1: Input to each iteration.
  % - *2: Output from each iteration, used as input to the next iteration.
  % - *3: Output totally (filled in during each iteration).
  % - n: iteration index. 1,2,...,N (will most likely stop before N).
  Bins1 = [dec2bin(Ons-1,N); dec2bin(Dcs-1,N)]; % width = N
  Covs1 = [false(nOn,1); true(nDc,1)];      % width = 1
  Nums1 = [Ons'-1; Dcs'-1];                 % width = 2^(n-1) in each iteration
  nVars1 = size(Bins1,1);
  twoexp = 2.^(N-1:-1:0)';                  % coeff to calculate Dons, [8;4;2;1]

  MINC = 100; % matrix increase size: Larger => faster for big tt, Smaller => slower for very big tt's.
  BinsNew = char('X'*ones(MINC,N));         % To prealloc an "empty" matrix (use 'X' for debug purpose)
  Bins2 = BinsNew;
  Covs2 = false(MINC,1);
  Nums2 = int16(-ones(MINC,1));
  % nVars2 initiated in each iteration.
  Bins3 = BinsNew;
  % Cov3 is not used, since we only move "uncovered" terms to *3.
  Nums3 = cell(MINC,1); % The content will differ in size - this saves memory.
  nVars3 = 0;
  for n=1:N % Main iteration loop. Will break in the middle (unless tt = ones(2^N,1))
    % If there is no terms to test, then we are done.
    if ~nVars1, break; end
    six = [find(~Covs1); find(Covs1)]; % Sort the variables so uncovered comes first
    Bins1 = Bins1(six,:);
    Covs1 = Covs1(six,:);
    Nums1 = Nums1(six,:);
    Dons1 = int16((Bins1 == '-')*twoexp);  % '-1-0' => [1 0 1 0]*[8;4;2;1] = 10
    Ones1 = uint16((Bins1 == '1')*twoexp); % '-1-0' => [0 1 0 0]*[8;4;2;1] = 4
    Nums2 = [Nums2,Nums2]; % We will double the number of atomic terms each time
    nVars2 = 0;
    
    % What do we have here, in each iteration?
    % n=1:
    %   *3 = empty
    %   Bins1 = [0100;1000;1010;1011;1100;1111;0111;1001;1110]
    %   Covs1 = [   0;   0;   0;   0;   0;   0;   1;   1;   1]
    %   Nums1 = [   4;   8;  10;  11;  12;  15;   7;   9;  14]
    %   Dons1 = [   0;   0;   0;   0;   0;   0;   0;   0;   0]
    %   Ones1 = [   4;   8;  10;  11;  12;  15;   7;   9;  14]
    % n=2:
    %   *3 = empty  1    2     3     4    5    6     7    8     9  10*   11*
    %   Bins1 = [-100;10-0; 101-; 1-11;1-00;100-; 1-10;10-1; 11-0;-111; 111-]
    %   Covs1 = [   0;   0;    0;    0;   1;   1;    1;   1;    1;   1;    1]
    %   Nums1 = [4,12;8,10;10,11;11,15;8,12; 8,9;10,14;9,11;12,14;7,15;14,15]
    %   Dons1 = [   8;   2;    1;    4;   4;   1;    4;   2;    2;   8;    1]
    %   Ones1 = [   4;   8;   10;   11;   8;   8;   10;   9;   12;   7;   14]
    % n=3:
    %   Bins3 = [  -100]
    %   Nums3 = {[4,12]} 1         2           3          4           5          6
    %   Bins1 = [     10--;     10--;       1-1-;      1--0;       1-1-;      1--0]
    %   Covs1 = [        0;        0;          0;         1;          1;         1]
    %   Nums1 = [8,9,10,11;8,9,10,11;10,11,14,15;8,10,12,14;10,11,14,15;8,10,12,14]
    %   Dons1 = [        3;        3;          5;         6:          5;         6]
    %   Dons1 = [        8;        8;         10;         8:         10;         8]
    % n=4 (before breaking due to ~nVars1):
    %   *1 = empty
    %   Bins3 = [  -100;       10--;         1-1-]
    %   Nums3 = {[4,12];[8,9,10,11];[10,11,14,15]}
    
    for ix0 = 1:nVars1-1
      % All *10 is current row from *1.
      Dons10 = Dons1(ix0);
      if Dons10 < 0, continue; end
      Bins10 = Bins1(ix0,:);
      Covs10 = Covs1(ix0);
      Nums10 = Nums1(ix0,:);
      Ones10 = Ones1(ix0);
      for ix1 = ix0+1:nVars1
        if Dons10 ~= Dons1(ix1), continue; end % Also works if Dons(ix1) = -1
        tmp = bitxor(Ones10,Ones1(ix1)); % binary version of Bins10 ~= Bins1(ix1,:)
        if tmp > 0 && bitand(tmp, tmp-1) > 0, continue; end % if at most 1 bit => continue.
        diffs = find(Bins10 ~= Bins1(ix1,:), 1); % find at most 1 pos. where they differs.
        if ~isempty(diffs) % Bins1(ix0,:) and Bins1(ix1,:) differs in only one pos
          nVars2 = nVars2 + 1;
          if length(Covs2) < nVars2 % grow *2 with MINC at a time
            Bins2 = [Bins2; BinsNew];
            Covs2 = [Covs2; false(MINC,1)];
            Nums2 = [Nums2; -ones(MINC,2^n)];
          end
          Bins2(nVars2,:) = Bins10;
          Bins2(nVars2,diffs) = '-';
          Covs2(nVars2) = false;%Covs10 & Covs1(ix1); % if both terms are covered, then this entire term is covered.
          % The above function (Covs10 & Covs1(ix1)) gave problem. Better to let Step2 handle the coverage.
          Nums2(nVars2,:) = [Nums10 Nums1(ix1,:)];
          Covs1([ix0 ix1]) = true;
          Covs10 = true;
        else % Those are equal => "remove" the first one
          Dons1(ix1) = -1;    % indicate a duplicate
          if ~Covs10, Covs1(ix1) = true; end % prohibit it from being moved to *3, if first one is
        end
      end
    end
    % Cover yet uncovered terms with *3: Uncovered *1 -> *3
    ixUC1 = find(~Covs1); % Number index for UnCovered terms
    if ~isempty(ixUC1)
      ixUCRng3 = nVars3 + (1:length(ixUC1)); % UnCovered, as index range for *3
      nVars3 = ixUCRng3(end);
      Bins3(ixUCRng3,:) = Bins1(ixUC1,:);
      Nums3(ixUCRng3) = num2cell(Nums1(ixUC1,:),2); % each row in Nums1 will be a cell element.
    end
    
    % Copy *2 -> *1: We don't want to test the already tested again, so only test the generated.
    Bins1 = Bins2(1:nVars2,:);
    Covs1 = Covs2(1:nVars2);
    Nums1 = Nums2(1:nVars2,:);
    nVars1 = nVars2;
  end
  
  if nVars1 > 0 % This can happend iff there are no '0' in tt.
    ixUCRng3 = nVars3 + (1:nVars1); % UnCovered, as index range for *3
    nVars3 = ixUCRng3(end);
    Bins3(ixUCRng3,:) = Bins1;
    Nums3(ixUCRng3) = num2cell(Nums1,2); % each row in Nums1 will be a cell element.
  end
  
  % Remove all no longer needed variables
  clear *1 *10 *2 ix* BinsNew Dcs MINC diffs n nDc nOn Ons six twoexp
  
  % Okay, what do we give to the next step? We give:
  % - Bins3: char-binary representation of the terms. E.g. [ -100 ;     10--   ;      1-1-    ]
  % - Nums3: cell array containing vectors            E.g. {[4,12]; [8,9,10,11]; [10,11,14,15]}
  %          with covered atomic terms.
  % - nVars3: integer counting terms.                 E.g. 3
  % - tt: original truth table (char).                E.g. '0000100-1-1110-1'
  % - N: number of boolean variables (tt = 2^N long)  E.g. 4

  %% Step 2: Reduce terms covered by other terms.
  
  % List which atomic terms are covered by which term
  A = false(nVars3,2^N);    % tt = '0000100-1-1110-1'
  for i = 1:nVars3          %  A = [0000100000001000
    A(i,Nums3{i}+1) = true; %       0000000011110000
  end                       %       0000000000110011] after the loop
  vcnt = sum(A,1);          % vcnt=[0000100011221011] vertical count of A.
  
  if verify
    % Check that no tt=='0' are covered, and all tt=='1' are covered.
    % This should never happend, independent of user input.
    % ("UNEXPECTED ERROR: ..." means there is an error in the script.)
    assert(all(vcnt(tt=='0') == 0), 'minTruthtable:zero_covered', 'UNEXPECTED ERROR: A ''0'' is covered.');
    assert(all(vcnt(tt=='1') >= 1), 'minTruthtable:ones_uncovered', 'UNEXPECTED ERROR: A ''1'' is not covered.');
  end
  
  % Remove from count and A all columns but thous with '1' in tt
  ttEq1 = tt=='1';        % A = [100010
  A = A(:,ttEq1);         %      011100
  vcnt = vcnt(:,ttEq1);   %      001101]
                          %vcnt=[112211]
  
  % Now we have a delicate problem:
  % cost = N-sum(Bins3=='-',2); ish
  % Minimize cost(vars4)
  %   vars4: linear index vector
  % s.t.    all(sum(A(vars4,:),1)>0)
  % In words: Find the selection of rows from A that covers all columns at least
  %           ones, and that minimizes the cost.
  
  if doexact
    % Ok. We have such an optimization function in the bottom. It's slow, but
    % let's use it.
    cost = N+1 - sum(Bins3(1:nVars3,:)=='-',2); % '-11-' => tree inputs. two AND + one OR.
    cost(cost==2) = 1; % a one-input AND gate = no gate.
    % ix = binopt(A,cost,maxcost). maxcost must be given, but do not really matter
    vars4 = find(binopt(A,cost,inf));
    clear A vcnt ttEq1 cost nVars3
  elseif doquick
    % ...or use the quicker heuristic function, also provided below.
    cost = N+1 - sum(Bins3(1:nVars3,:)=='-',2); % '-11-' => tree inputs. two AND + one OR.
    cost(cost==2) = 1; % a one-input AND gate = no gate.
    vars4 = find(quickopt(A,cost));
    clear A vcnt ttEq1 cost nVars3
  elseif ~isnan(dolevel)
    cost = N+1 - sum(Bins3(1:nVars3,:)=='-',2); % '-11-' => tree inputs. two AND + one OR.
    cost(cost==2) = 1; % a one-input AND gate = no gate.
    vars4 = find(level_opt(A,cost,inf,dolevel));
    clear A vcnt ttEq1 cost nVars3
  else
    
    % Detect all "Essential" prime implicant (that are alone on an atomic term).
    vars3 = (1:nVars3)'; % pointers of all *3 variables, used to index rows in A.
    vars4 = zeros(0,1);  % pointers of all *4 variables (points to rows in *3).
    vcntEq1 = find(vcnt==1);
    for c = vcntEq1(end:-1:1);  % look from right, to not "destroy" the index.
      if vcnt(c) == 1
        r = find(A(vars3,c),1); % which row (in vars3) contains the "true"?
        % r points at an essential term -> move this to *4
        vcnt(A(vars3(r),:)) = 0;        % Remove column by setting count to 0.
        vars4 = [vars4; vars3(r)];      % Add row to *4
        vars3 = vars3([1:r-1,r+1:end]); % Remove row from *3
      end
    end
    
    % If there happens to be any non-essential terms (rows) left, just pick the biggest first.
    % This cannot guarantee the best solution.
    cntDC = sum(Bins3(1:nVars3,:) == '-',2); % Number of '-' for each prod term.
    while any(vcnt > 0)
      % Remove the "removed" columns for real. Keep the row index for simplicity.
      A = A(:,vcnt > 0);
      vcnt = vcnt(vcnt > 0);
      hcnt = sum(A(vars3,:),2);       % Count used number of ones (indexing vars3)
      % Remove all terms that are by now covered by others.
      vars3 = vars3(hcnt > 0);        % Keep rows that can contribute.
      hcnt = hcnt(hcnt > 0);          % Update...
      %ccnt = A(vars3,:)*vcnt';        % Cover Count. ccnt(3) = how many '1' in A that will be removed if vars3(3) is selected.
      cntDC2 = cntDC(vars3);          % Count the number of '-' (indexing vars3)
      % Using count2 <=> smaller AND gates. Using count3 <=> fewer AND gates <=> smaller OR gate.
      [~,r] = max(hcnt*N+cntDC2);     % r = the vars3 index corresponding to biggest term, in some sense
      vcnt(A(vars3(r),:)) = 0;        % Remove columns by setting count to 0.
      vars4 = [vars4; vars3(r)];      % Add row to *4.
      vars3 = vars3([1:r-1,r+1:end]); % Remove row from *3.
    end
    vars4 = sort(vars4);
    clear A c vcnt count2 countEq1 i ixs nVars3 r ttEq1 vars3
  end
  
  % Okay, to the next step we leave:
  % - Bins3: (same as before).   E.g. [ -100 ;     10--   ;      1-1-    ]
  % - Nums3: (same as before).   E.g. {[4,12]; [8,9,10,11]; [10,11,14,15]}
  % - vars4: Index to the rows in Bins3. E.g. [1;2;3]
  
  %% Step 3: Finalize.
  Bins = Bins3(vars4,:);
  Nums = Nums3(vars4);
  inps = sum(Bins ~= '-', 2);
  inps = sum(inps .* (inps > 1),1) + length(inps)*(length(inps)>1);
  if nargout >= 4 || verify || nargout == 0 % calculate tt (out)
    ott = char('0'*ones(1,2^N));
    for i = 1:length(Nums)
      ott(Nums{i}+1) = '1';
    end
    if verify
      assert(all(ott(tt=='0') == '0'), 'minTruthtable:0_to_1', 'UNEXPECTED ERROR: A ''0'' become a ''1''.');
      assert(all(ott(tt=='1') == '1'), 'minTruthtable:1_to_0', 'UNEXPECTED ERROR: A ''1'' became a ''0''.');

      % Extra verification that Bins <=> Nums. Uncomment this if you suspect that the result is not correct.
      for i=1:size(Bins,1)
        Binsi = Bins(i,:);
        Zerosi = find(Binsi == '0');
        Onesi = find(Binsi == '1');
        Donsi = find(Binsi == '-');
        assert(length(Zerosi)+length(Onesi)+length(Donsi) == N, 'minTruthtable:CharInBins', 'UNEXPECTED ERROR: Unexpected character in Bins.');
        Numsi = sum(2.^(N-Onesi));
        for j=1:length(Donsi)
          Numsi = [Numsi (Numsi+2^(N-Donsi(j)))];
        end
        assert(length(Numsi) == length(Nums{i}), 'minTruthtable:LenNums_vs_Dons', 'UNEXPECTED ERROR: Length(Nums{%d}) ~= 2^(#-=%d).', i, length(Donsi));
        assert(all(sort(Numsi) == sort(Nums{i})), 'minTruthtable:Nums_vs_Dons', 'UNEXPECTED ERROR: Nums{%d} ~= Bins(%d,:).', i, i);
      end
    end
  end
  
  %% If no output required: Print result
  if nargout == 0
    if N==3 || N==4
      mat = [0 1 3 2; 4 5 7 6; 12 13 15 14; 8 9 11 10];
      head='/||\';
      if N==3, head='/\'; end
      fprintf('| Kmap (index to left, -:unused don''t care, =:used don''t care):\n')
      for y=1:2^(N-2)
        fprintf('|  %c %2d %2d %2d %2d %c ', head(y), mat(y,:), head(end+1-y));
        kmaptt = ott;
        kmaptt(tt=='0') = '.';
        kmaptt(tt=='-' & ott=='0') = '-';
        kmaptt(tt=='-' & ott=='1') = '=';
        fprintf(' %c %c %c %c %c %c\n', head(y), kmaptt(1+mat(y,:)), head(end+1-y));
      end
    end

    fprintf('| All terms:\n');
    for i=1:size(Bins,1)
      fprintf('| * T(%2d): "%s" <-> {', i,Bins(i,:))
      Tnums = Nums{i};
      if ~isempty(Tnums)
        fprintf('%d', Tnums(1));
        fprintf(' %d', Tnums(2:end));
      end
      fprintf('}\n');
    end
    if verify
      fprintf('| Relation "<->" has been verified\n');
    end
    fprintf('| y = ');
    ors = ''; % or-string, will change to ' + '
    for i=1:size(Bins,1)
      fprintf('%s',ors);
      ors = ' + ';
      ands = ''; % and-string, will change to '*'
      for j=1:size(Bins,2)
        x = size(Bins,2) - j;
        if Bins(i,j) == '0', fprintf('%s~x(%d)',ands,x); ands = '*';
        elseif Bins(i,j) == '1', fprintf('%sx(%d)',ands,x); ands = '*';
        end
      end
      if all(Bins(i,:) == '-'), fprintf('1'); end
    end
    if size(Bins,1) == 0, fprintf('0'); end
    fprintf(';\n| \n');
    fprintf('| Logical complexity: %d inputs\n', inps);
    fprintf('| \n');
    fprintf('| Input tt:  ''%s''\n', tt);
    fprintf('| Output tt: ''%s''\n', ott);
    clear Bins Nums inps ott
  end
end

%% TODOs:
% * Rewrite so the function can take several truth tables (a truth matrix), and co-optimize it.
%   - Suggested way: Run Step1 for each tt, and optimization shared logic in step 2.
% * Major change: Rewrite so 0 <= tt(i) <= 1, where tt(i) is double.
%   Then minimize cost = inps + sum(abs(tt - ott) / abs(1-ott-tt)).
%   In this way, the user can say "try to make this to a '1', but don't wast too many inputs on it".

function BIST()
  %% built-in self test
  digits = '0-1';
  fprintf('/== minTruthtable::BIST ==\n');
  for M=3:8 % number of variables
    ttn = ones(1,2^M); % ttn = numbered version of tt. 1 = first element in "digits"
    if M <= 3
      cnt = 3^(2^M);
    elseif M <= 6
      cnt = -1000;
    elseif M <= 7
      cnt = -100;
    else
      cnt = -9;
    end
    if cnt < 0
      fprintf('| /== %d inputs (vec len=2^%d=%d, #tests=%d out of 3^%d=%d) ==\n', M,M,2^M,-cnt,2^M,3^(2^M));
    else
      fprintf('| /== %d inputs (vec len=2^%d=%d, #tests = all 3^%d=%d) ==\n', M,M,2^M,2^M,cnt);
    end
    flags = {'e','12','8','5','3','2','1','q','Q',''};
    if M>=6, flags{1} = '14'; end
    sums = zeros(size(flags)); % sum(cost-cost_exact). Except sums(1) = sum(cost_exact)
    maxs = zeros(size(flags)); % max(cost-cost_exact)
    cntdiffs = zeros(size(flags));% count(cost>cost_exact)
    times = zeros(size(flags)); % sum(toc)
    fprintf('| | Progress (done=9): 0');
    next_i = round(abs(cnt)/9);
    for i=1:abs(cnt)
      if cnt<0, ttn = ceil(3*rand(1,2^M)); end
      tt = digits(ttn);
      t=tic; [~, inpe]=minTruthtable(tt,flags{1}); times(1) = times(1)+toc(t); sums(1)=sums(1)+inpe;
      for j=2:numel(flags)
        t=tic; [~, inp]=minTruthtable(tt,flags{j}); times(j)=times(j)+toc(t);
        sums(j) = sums(j) + (inp-inpe);
        maxs(j) = max(maxs(j), (inp-inpe));
        cntdiffs(j) = cntdiffs(j) + (inp>inpe);
        if strcmp(flags{1},'e')
          assert(inpe<=inp, 'minTruthtable:BIST:UNOPT', 'inpe=%d should be < inp(%s)=%d.',inpe,flags{j},inp);
        end
      end
      if cnt>=0
        % update the tt vector: least significant digit to the left.
        p=find(ttn~=3,1); % find position of first non-max digit. This is where the carry propagation stops.
        if ~isempty(p)
          ttn(1:p-1) = 1; % reset less significan digits.
          ttn(p) = ttn(p)+1; % increment the first non-max digit
        end
        assert(any(ttn>1), 'UNEXPECTED');
      end
      % print progress indicator
      if i==next_i
        foo = round(9*i/abs(cnt)); % 1, 2, ..., 9=at end
        fprintf('%d',foo);
        next_i = round((foo+1)*abs(cnt)/9);
      end
    end
    fprintf('\n'); % end of the progress bar
    fprintf('| | flag +-- sum* --+- max* -+-#diff*-+-- time --+\n')
    fprintf('| |%+5s |%9d |%7d |%7d |%8.3f s| (* compared to this row)\n',...
      ['''' flags{1} ''''],sums(1),maxs(1),cntdiffs(1),times(1));
    for j=2:numel(flags)
      fprintf('| |%+5s |%+9d |%+7d |%+7d |%8.3f s|\n',...
        ['''' flags{j} ''''],sums(j),maxs(j),cntdiffs(j),times(j));
    end
    fprintf('| \\--\n');
  end
  fprintf('\\--\n');
end

function ix = binopt(A, cost, maxcost)
  % Find those rows that _has_ to be included
  vcnt = sum(A,1); % count how many rows covers each column
  %assert(all(vcnt > 0), 'INTOPT:UNSOLVABLE');
  ix = any(A(:,vcnt==1),2); % logical index to rows in A
  
  % Find all that must be left. Exit if nothing is left
  % Keep lists of rows/cols that remains (left to use)
  cuse = find(~any(A(ix,:),1)); % Columns not covered by the required rows.
  if isempty(cuse), return, end
  ruse = find(~ix); % Each row in ix _has_ to be included. List the others.
  maxcost = maxcost - sum(cost(ix));
  if maxcost<0, ix=true(size(cost)); return; end
  %if maxcost<0, return; end
  
  % Remove unused columns
  hcnt = sum(A(ruse,cuse),2); % row index within ruse.
  ruse = ruse(hcnt > 0); % remove unused rows
  hcnt = hcnt(hcnt > 0);
  % Individually "optimize" each isolated column:
  isolated_cuse = find(~any(A(ruse(hcnt>1),cuse),1));
  if ~isempty(isolated_cuse)
    for i = isolated_cuse
      tmpix = find(A(ruse,cuse(i))); % index to ruse
      [~,minix] = min(cost(ruse(tmpix))); % find cheapest. index within tmpix.
      ix(ruse(tmpix(minix))) = true; % add it to result
      ruse(tmpix) = []; % remove rows
      hcnt(tmpix) = [];
    end
    cuse(isolated_cuse) = []; % remove the columns
  end
  if isempty(cuse), return, end
  if isempty(ruse), return, end
  % Todo: Identify isolated groups of column, that can be individually optimized.
  
  % Now, find the row that covers most, and "test" it on and off,
  hcnt = sum(A(ruse,cuse),2); % index within ruse
  [~,sel_row] = max(hcnt); % row index within ruse.
  sel_row2 = ruse(sel_row); % row index within A.
  covered_cols = A(sel_row2,cuse); % index to cuse, which are covered.
  ruse(sel_row) = [];      % remove the selected row.
  
  % test to have it on:
  cuse2 = cuse(~covered_cols);
  cost_on = cost(sel_row2); % to be updated
  if cost_on > maxcost % can't efford it
    ix_on = true(size(cost));
    cost_on = inf;
  elseif isempty(cuse2) % we are done
    ix_on = false(size(cost));
  else
    ix_on = ruse(binopt(A(ruse,cuse2),cost(ruse), maxcost-cost_on)); % index rows in A
    cost_on = cost_on + sum(cost(ix_on));
  end
  % update maxcost, so the "off" version does not have to find too expensive solutions.
  if maxcost > cost_on;
    maxcost = cost_on;
  end
  % test to have it off:
  ix_off = ruse(binopt(A(ruse,cuse),cost(ruse),maxcost)); % index rows in A
  cost_off = sum(cost(ix_off));
  
  % pick the best one:
  if cost_off < cost_on
    ix(ix_off) = true;
  else
    ix(ix_on) = true;
    ix(sel_row2) = true;
  end
end

function ix = level_opt(A, cost, maxcost, maxdepth)
  % Find those rows that _has_ to be included
  vcnt = sum(A,1); % count how many rows covers each column
  %assert(all(vcnt > 0), 'INTOPT:UNSOLVABLE');
  ix = any(A(:,vcnt==1),2); % logical index to rows in A
  
  % Find all that must be left. Exit if nothing is left
  % Keep lists of rows/cols that remains (left to use)
  cuse = find(~any(A(ix,:),1)); % Columns not covered by the required rows.
  if isempty(cuse), return, end
  ruse = find(~ix); % Each row in ix _has_ to be included. List the others.
  maxcost = maxcost - sum(cost(ix));
  if maxcost<0, ix=true(size(cost)); return; end
  %if maxcost<0, return; end
  
  if maxdepth<0
    % From here, skip the binary optimization, and use the quick-search instead.
    ix(ruse(quickopt(A(ruse,cuse),cost(ruse)))) = true;
    return
  end
  
  % Remove unused columns
  hcnt = sum(A(ruse,cuse),2); % row index within ruse.
  ruse = ruse(hcnt > 0); % remove unused rows
  hcnt = hcnt(hcnt > 0);
  % Individually "optimize" each isolated column:
  isolated_cuse = find(~any(A(ruse(hcnt>1),cuse),1));
  if ~isempty(isolated_cuse)
    for i = isolated_cuse
      tmpix = find(A(ruse,cuse(i))); % index to ruse
      [~,minix] = min(cost(ruse(tmpix))); % find cheapest. index within tmpix.
      ix(ruse(tmpix(minix))) = true; % add it to result
      ruse(tmpix) = []; % remove rows
      hcnt(tmpix) = [];
    end
    cuse(isolated_cuse) = []; % remove the columns
  end
  if isempty(cuse), return, end
  if isempty(ruse), return, end
  % Todo: Identify isolated groups of column, that can be individually optimized.
  
  % Now, find the row that covers most, and "test" it on and off,
  hcnt = sum(A(ruse,cuse),2); % index within ruse
  [~,sel_row] = max(hcnt); % row index within ruse.
  sel_row2 = ruse(sel_row); % row index within A.
  covered_cols = A(sel_row2,cuse); % index to cuse, which are covered.
  ruse(sel_row) = [];      % remove the selected row.
  
  % test to have it on:
  cuse2 = cuse(~covered_cols);
  cost_on = cost(sel_row2); % to be updated
  if cost_on > maxcost % can't efford it
    ix_on = true(size(cost));
    cost_on = inf;
  elseif isempty(cuse2) % we are done
    ix_on = false(size(cost));
  else
    ix_on = ruse(level_opt(A(ruse,cuse2),cost(ruse), maxcost-cost_on, maxdepth-1)); % index rows in A
    cost_on = cost_on + sum(cost(ix_on));
  end
  % update maxcost, so the "off" version does not have to find too expensive solutions.
  if maxcost > cost_on;
    maxcost = cost_on;
  end
  % test to have it off:
  ix_off = ruse(level_opt(A(ruse,cuse),cost(ruse),maxcost,maxdepth-1)); % index rows in A
  cost_off = sum(cost(ix_off));
  
  % pick the best one:
  if cost_off < cost_on
    ix(ix_off) = true;
  else
    ix(ix_on) = true;
    ix(sel_row2) = true;
  end
end

function ix = quickopt(A, cost)
  % Recursion termination, v1:
  if isempty(A)
    ix = false(size(A,1),1);
    return
  end
  % Find those rows that _has_ to be included
  vcnt = sum(A,1); % count how many rows covers each column
  assert(all(vcnt > 0), 'INTOPT:UNSOLVABLE');
  ix = any(A(:,vcnt==1),2); % logical index to rows in A
  
  % Remove the corresponding rows and columns (by adding index of remaining)
  ruse = find(~ix); % Each row in ix _has_ to be included. List the others.
  cuse = find(~any(A(ix,:),1)); % Columns not covered by the required rows.
  
  % Recursion termination, v2:
  if isempty(cuse)
    return
  end
  
  % Now, try to find a suitable cost model, and then take the minimum of this.
  % * column weight = 1./vcnt.
  % * row weight = sum(weight of covered columns).
  % * find min(cost ./ row weight)
  row_weight = A(ruse,cuse) * (1./vcnt(cuse)'); % matrix mul = sum(1./vcnt(selected))
  cost2 = cost(ruse) ./ row_weight; % cost per weighted cover... ish...
  [~,sel_row] = min(cost2); % indexing ruse
  
  % Remove this row and its corresponding columns
  cuse(A(ruse(sel_row), cuse)) = []; % remove columns
  ix(ruse(sel_row)) = true; % save the selected to ix
  ruse(sel_row) = [];       % remove
  
  % recurse
  % TODO: loopify, to get rid of the recursion
  ix(ruse(quickopt(A(ruse,cuse), cost(ruse)))) = true;
end
