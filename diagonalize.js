function norm(vector){
	var sum = 0;

	for(var i = 0; i < vector.length; i++){
		sum += math.re(math.multiply(vector[i], math.conj(vector[i])));
	}

	return math.sqrt(sum);
}

function norm_squared(vector){
	var sum = 0;

	for(var i = 0; i < vector.length; i++){
		sum += math.re(math.multiply(vector[i], math.conj(vector[i])));
	}

	return sum;
}

function givens_rotate(vector, c, s, x, y){
	var x_1 = vector[x]*c - vector[y]*s;
	var y_1 = vector[x]*s + vector[y]*c;

	vector[x] = x_1;
	vector[y] = y_1;

	return vector;
}


function givens_rotation(n, c, s, x, y){
	var output = math.identity(n).toArray();

	output[x][x] = c;
	output[y][y] = c;
	output[x][y] = -s;
	output[y][x] = s;

	return output;
}

function hessenberg(matrix){
	if(matrix.length <= 2){
		return {unitary: math.identity(matrix.length), result: matrix};
	}

	var m1 = math.subset(matrix, math.index(math.range(1, matrix.length), math.range(0, matrix.length)));
	var first_column = math.subset(matrix, math.index(math.range(1, matrix.length), math.range(0, 1)));
	var first_vector = [[norm(first_column)]];

	for(var i = 1; i < matrix.length - 1; i++){
		first_vector[i] = [0];
	}
	
	if(Math.abs(math.re(first_column[0])) < math.config.epsilon && Math.abs(math.im(first_column[0])) < math.config.epsilon){
		var w = math.subtract(first_vector, first_column);
	} else {
		var w = math.add(first_vector, math.multiply(first_column[0][0], 1.0/math.abs(first_column[0][0]), first_column));
	}

	var w_adj = math.ctranspose(w);
	var V = math.subtract(math.identity(matrix.length - 1), math.multiply(2, 1.0/norm_squared(w_adj), w, w_adj));
	var U = math.identity(matrix.length);
	U = math.subset(U, math.index(math.range(1, matrix.length), math.range(1, matrix.length)), V);
	var next_matrix = math.multiply(U, matrix, U);

	var next_result = hessenberg(math.subset(next_matrix, math.index(math.range(1, matrix.length), math.range(1, matrix.length))).toArray());

	var new_U = math.identity(matrix.length);
	new_U = math.subset(new_U, math.index(math.range(1, matrix.length), math.range(1, matrix.length)), next_result.unitary);

	var next_matrix = math.subset(next_matrix, math.index(math.range(1, matrix.length), math.range(1, matrix.length)), next_result.result);

	return {unitary: math.multiply(U, new_U), result: next_matrix};
}

function tridiag_to_matrix(diag, subdiag){
	var output = math.identity(diag.length).toArray();

	for(var i = 0; i < diag.length - 1; i++){
		output[i][i] = diag[i];
		output[i + 1][i] = subdiag[i];
		output[i][i + 1] = math.conj(subdiag[i]);
	}

	output[diag.length - 1][diag.length - 1] = diag[diag.length - 1];

	return output;
}

//Not correct
//I must *conjugate* by the givens rotation below in order to apply the implicit q theorem
function implicit_qr_step(diag, subdiag){
	var x = diag[0];
	var y = subdiag[0];

	var d = math.sqrt(x*x + y*y);
	var c = x/d;
	var s = -y/d;

	diag[0] = d;
	subdiag[0] = 0;
	diag[1] = y*s + c*diag[1];
	var bulge = -s*subdiag[1];
	subdiag[1] *= c;

	return {diag: diag, subdiag: subdiag, bulge: bulge};
}
