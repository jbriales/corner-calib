% Generated using GBSolver generator Copyright Martin Bujnak,
% Zuzana Kukelova, Tomas Pajdla CTU Prague 2008.
% 
% Please refer to the following paper, when using this code :
%     Kukelova Z., Bujnak M., Pajdla T., Automatic Generator of Minimal Problem Solvers,
%     ECCV 2008, Marseille, France, October 12-18, 2008

function [x y z] = sw_Zhou_alt(d1_1, d2_1, d3_1, b_1)

	% precalculate polynomial equations coefficients
	c(1) = d1_1(1);
	c(2) = d1_1(2);
	c(3) = d1_1(3);
	c(4) = -1;
	c(5) = d2_1(1);
	c(6) = d2_1(2);
	c(7) = d2_1(3);
	c(8) = -1;
	c(9) = d3_1(1);
	c(10) = 2*b_1(3);
	c(11) = d3_1(2);
	c(12) = 2*b_1(2);
	c(13) = 2*b_1(1);
	c(14) = d3_1(3);
	c(15) = 0;

	M = zeros(26, 34);
	ci = [17, 42, 119, 144, 221, 370, 395, 472, 625];
	M(ci) = c(1);

	ci = [69, 94, 171, 196, 273, 422, 447, 524, 677];
	M(ci) = c(2);

	ci = [251, 276, 301, 326, 351, 552, 577, 602, 755];
	M(ci) = c(3);

	ci = [667, 692, 717, 742, 767, 786, 811, 836, 859];
	M(ci) = c(4);

	ci = [22, 47, 124, 149, 226, 373, 398, 475, 626];
	M(ci) = c(5);

	ci = [74, 99, 176, 201, 278, 425, 450, 527, 678];
	M(ci) = c(6);

	ci = [256, 281, 306, 331, 356, 555, 580, 605, 756];
	M(ci) = c(7);

	ci = [672, 697, 722, 747, 772, 789, 814, 839, 860];
	M(ci) = c(8);

	ci = [26, 129, 154, 231, 376, 401, 478, 627];
	M(ci) = c(9);

	ci = [52, 155, 180, 257, 402, 427, 504, 653];
	M(ci) = c(10);

	ci = [78, 181, 206, 283, 428, 453, 530, 679];
	M(ci) = c(11);

	ci = [156, 233, 258, 309, 480, 505, 556, 705];
	M(ci) = c(12);

	ci = [182, 259, 284, 335, 506, 531, 582, 731];
	M(ci) = c(13);

	ci = [260, 311, 336, 361, 558, 583, 608, 757];
	M(ci) = c(14);

	ci = [676, 727, 752, 777, 792, 817, 842, 861];
	M(ci) = c(15);


	Mr = rref(M);  % replace me with a MEX

	A = zeros(8);
	amcols = [34 33 32 31 30 29 28 24];
	A(1, 4) = 1;
	A(2, 7) = 1;
	A(3, :) = -Mr(25, amcols);
	A(4, :) = -Mr(24, amcols);
	A(5, :) = -Mr(22, amcols);
	A(6, :) = -Mr(20, amcols);
	A(7, :) = -Mr(19, amcols);
	A(8, :) = -Mr(12, amcols);

	[V D] = eig(A);
	sol =  V([4, 3, 2],:)./(ones(3, 1)*V(1,:));

	if (find(isnan(sol(:))) > 0)
		
		x = [];
		y = [];
		z = [];
	else
		
		I = find(not(imag( sol(1,:) )));
		x = sol(1,I);
		y = sol(2,I);
		z = sol(3,I);
	end
end
