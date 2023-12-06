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
