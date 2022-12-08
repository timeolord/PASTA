### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ f1d327bf-42ce-480e-b98c-42ac5ba67235
begin
	using Chain
	using Dictionaries
	using Random
	using StatsBase
	using BenchmarkTools
	using Unzip
	using Parameters
	using ProfileSVG
	using FASTX
	using BioSequences
	using Kmers
	using BioAlignments
	using Lazy
end

# ╔═╡ 4a7dc282-7355-4427-ad0a-efc1c3f9faad
begin 
	sequence = @chain read("seq.fa", String) begin
		LongDNA{2}(_)
	end
end

# ╔═╡ 74f5561d-2039-40f0-88a2-b19fe555c3a2
function to_nucleotide(character::T) where T<:AbstractChar
	if character == 'A'
		Nucleotide(1)
	elseif character == 'T'
		Nucleotide(2)
	elseif character == 'G'
		Nucleotide(3)
	else character == 'C'
		Nucleotide(4)
	end
end

# ╔═╡ 9e7f6e3a-b646-4fe4-bf8b-26159f816505
begin
	probabilites = @chain read("probs.txt", String) begin
		split(" ")
		@aside pop!(_)
		parse.(Float64, _)
	end
end

# ╔═╡ 4af62ee8-f649-400f-8182-1ac42ded6ae5
begin
	const nucleotides = [DNA_A, DNA_T, DNA_G, DNA_C]
	#EveryKmer{DNAKmer{hp.word_size}}(sequence) |> collect
end

# ╔═╡ 4d4578e4-4c67-45bd-a360-70193aae5d5e
prob_seq = (sequence = sequence, probabilites = probabilites)

# ╔═╡ 5625822c-2030-4303-88d1-0c5e95cae278
function windows(sequence, step_size, window_size)
    ((@view sequence[i:i+window_size-1]) for i in 1:step_size:length(sequence)-window_size+1)
end

# ╔═╡ 4c4979f3-dd90-4744-8c55-8eca1f9ccbea
# ╠═╡ disabled = true
#=╠═╡
struct ProbabilisticNucleotide
	nucleotide::DNA
	probability::Float64
end
  ╠═╡ =#

# ╔═╡ 6db80ce9-4c47-4672-b0c6-7b81424ceb71
#=╠═╡
Base.show(io::IO, x::ProbabilisticNucleotide) = print(io, x.nucleotide, " ", x.probability)
  ╠═╡ =#

# ╔═╡ 57c5982a-1df9-4cc2-9d46-537cd31e2d68
#=╠═╡
function Base.show(io::IO, ::MIME"text/plain", x::ProbabilisticNucleotide)
	print(io, x.nucleotide, " ", x.probability)
end
  ╠═╡ =#

# ╔═╡ f4857a38-02bb-4c02-89d4-aa4906229b3e
function to_int(DNA) 
	if DNA == DNA_A
		1
	elseif DNA == DNA_T
		2
	elseif DNA == DNA_G
		3
	else
		4
	end
end

# ╔═╡ 9edfc292-0aa5-4e07-8381-d7d02ec1a09d
function probability(nucleotide, probability)
	prob_rest = (1 - probability)/3
	result = [prob_rest for _ in nucleotides]
	result[to_int(nucleotide)] = probability
	result
end

# ╔═╡ 62d5b131-fd26-4225-8420-cd8b42a7200b
function sample_range(sequence, probabilities)
	map((character, prob)->
		sample(nucleotides, ProbabilityWeights(probability(character, prob))), sequence, probabilities)
end

# ╔═╡ eb97c760-a32f-42d5-8c7b-f9138b576fe2
# ╠═╡ disabled = true
#=╠═╡
function permutations(hyperparameters)
	[[Nucleotide(a) for a in Tuple(index)] for index in CartesianIndices(tuple(fill(4, hyperparameters.word_size)...))]
end
  ╠═╡ =#

# ╔═╡ 61b0eb60-f1c4-4067-9d70-ee181dd64a36
# ╠═╡ disabled = true
#=╠═╡
@with_kw struct Hyperparameters
	query_length_range::UnitRange{Int} = 10:100
	word_size::Int = 5
	permutation_threshold::Float64 = 0.9
	score_threshold::Int = 18
	data_base_score::AffineGapScoreModel{Int64} =  AffineGapScoreModel(EDNAFULL, gap_open=-500, gap_extend=-100)
	scoring_matrix = EDNAFULL
end
  ╠═╡ =#

# ╔═╡ 3220651e-c49b-4cc6-b76d-26790451affc
# ╠═╡ disabled = true
#=╠═╡
function prob_word_score2(query, word, probability, hyperparameters)
	sum(hyperparameters.scoring_matrix[q, w] * p for (q, w, p) in zip(query, word, probability))
end
  ╠═╡ =#

# ╔═╡ b5d3d606-91d6-454c-846b-81c744914042
# ╠═╡ disabled = true
#=╠═╡
function find_high_scoring_words2(query)
	filter(x->non_prob_word_score(query, x) > score_threshold, keys(blast_db))
end
  ╠═╡ =#

# ╔═╡ 699cfd94-5271-4aaa-8bae-5772324f58a3
# ╠═╡ disabled = true
#=╠═╡
struct SeedIterator
	length::Int
	query::DNAKmer{word_size}
end
  ╠═╡ =#

# ╔═╡ 3559f2e9-9a00-44f7-bc7e-1c05eab7af27
#=╠═╡
function Base.iterate(x::SeedIterator)
	index = 0
	range = blast_db[x.query][index]
	seq = @view sequence[range]
	prob = @view probabilites[range]
	while prob_word_score(x.query, seq, prob) < score_threshold
		index += 1
		if index > x.length
			return nothing
		end
		range = blast_db[x.query][index]
		seq = @view sequence[range]
		prob = @view probabilites[range]
	end
	(range, index + 1)
end
  ╠═╡ =#

# ╔═╡ ef8a4d03-e6fb-4e8e-9a8c-54f0ac2d8c7f
#=╠═╡
function Base.iterate(x::SeedIterator, i)
	index = i
	range = blast_db[x.query][index]
	seq = @view sequence[range]
	prob = @view probabilites[range]
	while prob_word_score(x.query, seq, prob) < score_threshold
		index += 1
		if index > x.length
			return nothing
		end
		range = blast_db[x.query][index]
		seq = @view sequence[range]
		prob = @view probabilites[range]
	end
	(range, index + 1)
end
  ╠═╡ =#

# ╔═╡ fa9cfa41-26d6-4817-a063-ca4e46b91f44
# ╠═╡ disabled = true
#=╠═╡
function find_seeds2(query)
		SeedIterator(length(blast_db[query]), query)
end
	
  ╠═╡ =#

# ╔═╡ fdff2519-1b53-45a4-ad00-8bb4a382bf9d
# ╠═╡ disabled = true
#=╠═╡
function find_seeds3(query)
		function is_seed(range)
			seq = @view sequence[range]
			prob = @view probabilites[range]
			if prob_word_score(query, seq, prob) > score_threshold
				true
			else
				false
			end
		end
		filter(is_seed, @lazy blast_db[query])
	end
  ╠═╡ =#

# ╔═╡ 58be55eb-d107-4d20-b2f2-a084f83ae338
# ╠═╡ disabled = true
#=╠═╡
function prob_word_score2(query, word, probability)
	sum(Int(a == b) for (a, b) in zip(query.data, word.data))
end
  ╠═╡ =#

# ╔═╡ 9ecd59ea-174b-4ad0-b943-d05566088bc3
#@btime PBLAST(query) samples=10

# ╔═╡ f7470de9-4a5a-4395-b2c5-ccf91734824b
#range

# ╔═╡ 623e3c50-8d86-418f-98c2-7a9c9da769bd
#pairalign(LocalAlignment(), query, subject, score_model)

# ╔═╡ 0035cbb3-6f3b-4e04-aa2c-d125b8baffbe
#@benchmark find_alignment($query, $subject, $probs) samples=5000

# ╔═╡ 5153c32d-2881-4461-9077-3476ced37798
# ╠═╡ disabled = true
#=╠═╡
@benchmark prob_word_score2($query, $subject, $probs) samples=5000
  ╠═╡ =#

# ╔═╡ ea42ca53-d31f-4370-a350-1f1f9fde8f41
# ╠═╡ disabled = true
#=╠═╡
function find_alignments2(query, word_index, seq_indices)
	filtered_indices = filter!(x->first(x) + word_index + length(query) < length(sequence), seq_indices)
	alignments = Array{Tuple{Float32, UnitRange{Int64}}}(undef, length(filtered_indices))
	Threads.@threads for i in 1:length(filtered_indices)
		seq_index = filtered_indices[i]
		alignments[i] = find_alignment(query, word_index, seq_index)
	end
	#alignments = map(seq_index -> find_alignment(query, word_index, seq_index), filtered_indices)
	maximum(alignments, init = (-Inf, []))
end
  ╠═╡ =#

# ╔═╡ 6e2ed658-bc97-4749-8b84-243c534c7d34
# ╠═╡ disabled = true
#=╠═╡
function find_alignments2(query, word_index, seq_indices)
	filtered_indices = filter(x->first(x) + word_index + length(query) < length(sequence), seq_indices)
	alignments = lazymap(seq_index -> find_alignment(query, word_index, seq_index), filtered_indices)
	maximum(alignments, init = (-Inf, []))
end
  ╠═╡ =#

# ╔═╡ 59a677a6-c2f4-4921-bcc8-dbe27242fc44
# ╠═╡ disabled = true
#=╠═╡
function PBLAST_max(query)
	words_list = to_words(query)
	seeds = map(
		((word_index, query_word),) -> 
		(word_index, find_seeds(query_word)), words_list)
	alignments = map(
		((word_index, seq_indices),) ->
		find_alignments(query, word_index, seq_indices), seeds)
	(a, b) = maximum(alignments)
	(score = a, range = b, sequence = sequence[b])
end
  ╠═╡ =#

# ╔═╡ ed31ab1c-70ea-4f6d-9824-f4cb6f8dc48d
# ╠═╡ disabled = true
#=╠═╡
function PBLAST_max2(query)
	words_list = to_words(query) |> collect
	alignments_list = Array{Tuple{Float64, UnitRange{Int64}}}(undef, length(words_list))
	Threads.@threads for i in 1:length(words_list)
		(word_index, query_word) = words_list[i]
		seeds = find_seeds(query_word)
		alignments_list[i] = find_alignments2(query, word_index, seeds)
	end
	(a, b) = maximum(alignments_list)
	(score = a, range = b, sequence = sequence[b])
end
  ╠═╡ =#

# ╔═╡ d13c8a85-e7d5-4a81-a70a-57b8ab6206a0
# ╠═╡ disabled = true
#=╠═╡
function PBLAST_max2(query)
	words_list = to_words(query)
	seeds = lazymap(
		((word_index, query_word),) -> 
		(word_index, find_seeds3(query_word)), words_list)
	alignments = lazymap(
		((word_index, seq_indices),) ->
		find_alignments2(query, word_index, seq_indices), seeds)
	(a, b) = maximum(alignments)
	(score = a, range = b, sequence = sequence[b])
end
  ╠═╡ =#

# ╔═╡ 0c6c200f-732d-4831-95f4-5f977e7a2b48
function perfect_match(query_range, result)
	query_range == result.range
end

# ╔═╡ 81f8a645-53ff-4e67-b6cf-25e3343c4b04
function start_match(query_range, result)
	first(query_range) == first(result.range)
end

# ╔═╡ e92a87b0-d37a-4b1b-9740-fe89a963cb16
function end_match(query_range, result)
	last(query_range) == last(result.range)
end

# ╔═╡ ac9c68b3-0a7a-468d-9b84-9c4a668c3c32
#Hyperparameters
begin
	const scoring_matrix::SubstitutionMatrix{DNA, Int64} = EDNAFULL
	#const scoring_matrix = DichotomousSubstitutionMatrix(1, -1)
	const gap_cost = -5
	const score_model::AffineGapScoreModel{Int64} =  AffineGapScoreModel(EDNAFULL, gap_open=gap_cost, gap_extend=gap_cost)
	const query_length_range::UnitRange{Int} = 100:500
	const word_size::Int = 5
	const permutation_threshold::Float64 = 0.9
	const T_score_threshold::Int = 18
	const S_score_cutoff::Int = 0
	const query_mutations_percent = 0
	const close_match_percent = 0.25
end

# ╔═╡ ccc968fa-f75f-4a30-90d1-cf0c97288ae0
function to_words(sequence)
	EveryKmer{DNAKmer{word_size}}(sequence)
end

# ╔═╡ c8948347-6cd9-41ab-aaf3-2e108ce02967
function generate_blast_db()
	blast_db = Dictionary{DNAKmer{word_size}, Vector{UnitRange{Int64}}}()
	for (index, window) in EveryKmer{DNAKmer{word_size}}(sequence)
		index_list = get!(blast_db, window, [])
		push!(index_list, index:index + word_size - 1)
	end
	blast_db
end

# ╔═╡ d43926d5-d305-4e12-bfc3-d72d6aa8c000
const blast_db = generate_blast_db()

# ╔═╡ c0158950-9a2e-42d3-8d82-3249cad660a2
function non_prob_word_score(query, word)
	sum(scoring_matrix[q, w] for (q, w) in zip(query, word))
end

# ╔═╡ c485200c-8637-45f9-85e6-df399758a09e
if T_score_threshold == 18 && word_size == 5
	function find_high_scoring_words(query)
		query
	end
else
	function find_high_scoring_words(query)
		filter(x->non_prob_word_score(query, x) > T_score_threshold, keys(blast_db))
	end
end

# ╔═╡ a1dc1cb0-58c5-44cf-8256-33d19c97b231
function prob_word_score(query, word, probability)
	function scoring(q, w, p)
		sum = 0
		for nucleotide in nucleotides
			if nucleotide == w
				sum += scoring_matrix[q, nucleotide] * p
			else
				sum += scoring_matrix[q, nucleotide] * (1 - p)/3
			end
		end
		sum
	end
	sum(scoring(q, w, p) for (q, w, p) in zip(query, word, probability))
end

# ╔═╡ 8b7aab0c-3690-4549-b840-e704ad1b70e4
if T_score_threshold == 18 && word_size == 5
	function find_seeds(query)
		seeds = []
		for range in blast_db[query]
			seq = sequence[range]
			prob = probabilites[range]
			if prob_word_score(query, seq, prob) > T_score_threshold
				push!(seeds, range)
			end
		end
		seeds
	end
else
	function find_seeds(query)
		high_score_words = find_high_scoring_words(query)
	
		blast_seqs = Iterators.flatten(blast_db[query] for query in high_score_words)
		@chain blast_seqs begin
			map(x->(sequence[x], probabilites[x], x), _)
			map(a->(prob_word_score(query, a[1], a[2]), last(a)), _)
			filter!(x->first(x) > T_score_threshold, _)
			last.()
		end
	end
end

# ╔═╡ 554206ee-5a38-47c3-9c68-d92fafbd6401
function find_alignment(query, seed_index, sequence_index)
	seq_start = first(sequence_index) + (seed_index - 1)
	seq_end = first(sequence_index) + (seed_index - 1) + length(query) - 1
	seq_range = seq_start:seq_end
	subject = @view sequence[seq_range]
	probs = @view probabilites[seq_range]
	(prob_word_score(query, subject, probs), seq_range)
	#alignment_ = pairalign(GlobalAlignment(), query, subject, score_model)
	#(BioAlignments.score(alignment_), seq_range)
end

# ╔═╡ 69559e82-65ee-4986-95f7-eead25f83e56
function find_alignments(query, word_index, seq_indices)
	filtered_indices = filter!(x->first(x) + word_index + length(query) < length(sequence), seq_indices)
	alignments = map(seq_index -> find_alignment(query, word_index, seq_index), filtered_indices)
	maximum(alignments, init = (-Inf, 1:1))
end

# ╔═╡ 4870eeeb-cd46-41d4-98e7-f87af65b3f21
function PBLAST(query; results_n=50)
	words_list = to_words(query) |> collect
	alignments_list = Array{Tuple{Float64, UnitRange{Int64}}}(undef, length(words_list))
	Threads.@threads for i in 1:length(words_list)
		(word_index, query_word) = words_list[i]
		seeds = find_seeds(query_word)
		alignments_list[i] = find_alignments(query, word_index, seeds)
	end
	top_list = Iterators.take(sort!(alignments_list, rev=true), results_n)
	map(((a, b),) -> (score = a, range = b, sequence = sequence[b]), top_list)
end

# ╔═╡ 88fe623f-6163-4e8e-8e34-5ce967c84ae2
function PBLAST_max(query)
	words_list = to_words(query) |> collect
	alignments_list = Array{Tuple{Float64, UnitRange{Int64}}}(undef, length(words_list))
	Threads.@threads for i in 1:length(words_list)
		(word_index, query_word) = words_list[i]
		seeds = find_seeds(query_word)
		alignments_list[i] = find_alignments(query, word_index, seeds)
	end
	(a, b) = maximum(alignments_list)
	(score = a, range = b, sequence = sequence[b])
end

# ╔═╡ 8fa9fd5c-ebd9-46c5-9695-471aef7b156c
function PBLAST_list(queries; results_n)
	if results_n == 1
		[PBLAST_max(query) for query in queries]
	else 
		[PBLAST(query, results_n = results_n) for query in queries]
	end
end

# ╔═╡ 0100a578-9f44-41a0-b0f0-dc2332d42325
# ╠═╡ disabled = true
#=╠═╡
function find_alignment2(query, seed_index, sequence_index)
	seq_start = first(sequence_index) + (seed_index - 1)
	seq_end = first(sequence_index) + (seed_index - 1) + length(query) - 1
	seq_range = seq_start:seq_end
	subject = @view sequence[seq_range]
	probs = @view probabilites[seq_range]
	#(prob_word_score(query, subject, probs), seq_range)
	alignment_ = pairalign(GlobalAlignment(), query, subject, score_model)
	(BioAlignments.score(alignment_), seq_range)
end
  ╠═╡ =#

# ╔═╡ 99b513ff-9928-4aba-8818-502c9136cbb0
#=╠═╡
function find_alignments2(query, word_index, seq_indices)
	filtered_indices = filter!(x->first(x) + word_index + length(query) < length(sequence), seq_indices)
	alignments = map(seq_index -> find_alignment2(query, word_index, seq_index), filtered_indices)
	maximum(alignments, init = (-Inf, []))
end
  ╠═╡ =#

# ╔═╡ 7e3d62bb-0892-4e5a-991c-0fbbfd7b2ae9
#=╠═╡
function PBLAST_max2(query)
	words_list = to_words(query) |> collect
	alignments_list = Array{Tuple{Float64, UnitRange{Int64}}}(undef, length(words_list))
	Threads.@threads for i in 1:length(words_list)
		(word_index, query_word) = words_list[i]
		seeds = find_seeds(query_word)
		alignments_list[i] = find_alignments2(query, word_index, seeds)
	end
	(a, b) = maximum(alignments_list)
	(score = a, range = b, sequence = sequence[b])
end
  ╠═╡ =#

# ╔═╡ 373d13c9-c11f-4a03-8ea5-5aef62471c4c
# ╠═╡ disabled = true
#=╠═╡
function find_alignment3(query, seed_index, sequence_index)
	seq_start = first(sequence_index) + (seed_index - 1)
	seq_end = first(sequence_index) + (seed_index - 1) + length(query) - 1
	seq_range = seq_start:seq_end
	subject = @view sequence[seq_range]
	probs = @view probabilites[seq_range]
	
	seed_start = seed_index
	seed_end = seed_index + word_size - 1
	#q = @view query[seed_start:seed_end]
	q = Array{Union{Nothing, DNA}}(nothing, length(query))
	#q::Array{Union{Nothing, DNA}} = fill(nothing, length(query))
	for index in seed_start:seed_end
		q[index] = query[index] 
	end
	score = prob_word_score(q, subject, probs)
	
	while score > S_score_cutoff
		if seed_start != 1
			seed_start -= 1
		end
		if seed_end != length(query)
			seed_end += 1
		end
		q = @view query[seed_start:seed_end]
		score = prob_word_score(q, subject, probs)
	end
	
	(score, seq_range)
end
  ╠═╡ =#

# ╔═╡ 22bca925-c356-4a4c-a963-d2ffb42af719
#=╠═╡
function find_alignments3(query, word_index, seq_indices)
	filtered_indices = filter!(x->first(x) + word_index + length(query) < length(sequence), seq_indices)
	alignments = map(seq_index -> find_alignment3(query, word_index, seq_index), filtered_indices)
	maximum(alignments, init = (-Inf, []))
end
  ╠═╡ =#

# ╔═╡ ae3725c3-0f0c-431c-8217-ca55c66aac05
#=╠═╡
function PBLAST_max3(query)
	words_list = to_words(query) |> collect
	alignments_list = Array{Tuple{Float64, UnitRange{Int64}}}(undef, length(words_list))
	Threads.@threads for i in 1:length(words_list)
		(word_index, query_word) = words_list[i]
		seeds = find_seeds(query_word)
		alignments_list[i] = find_alignments3(query, word_index, seeds)
	end
	(a, b) = maximum(alignments_list)
	(score = a, range = b, sequence = sequence[b])
end
  ╠═╡ =#

# ╔═╡ 4980de9a-b1c8-49b0-9db2-bc38f76b13ee
function close_match(query_range, result)
	len = length(query_range) * close_match_percent
	(first(result.range) in (first(query_range) - len - 1):(first(query_range) + len)) || (last(result.range) in (last(query_range) - len - 1):(last(query_range) + len))
end

# ╔═╡ cf34efa1-87ab-4094-a667-cf679f6edecf
function matches(queries, results)
	function match(query, result)
		if perfect_match(query, result)
			:P
		elseif start_match(query, result)
			:S
		elseif end_match(query, result)
			:E
		elseif close_match(query, result)
			:C
		else
			:N
		end
	end
	map(match, queries, results)
end

# ╔═╡ 5baedc9a-654b-4b96-8150-324ecf94c4ee
function mutate!(sequence)
	function insertion(sequence, index)
		insert!(sequence, index, rand(nucleotides))
	end
	function deletion(sequence, index)
		deleteat!(sequence, index)
	end
	function substitution(sequence, index)
		sequence[index] = rand(nucleotides)
	end
	
	mutation_functions = [insertion, deletion, substitution]
	
	mutations = rand(0:floor(length(sequence) * (query_mutations_percent)))
	if mutations != 1:1
		for _ in 1:mutations
			index = rand(1:length(sequence))
			rand(mutation_functions)(sequence, index)
		end
	end
	sequence
end

# ╔═╡ 8006ffbb-e654-4aff-9666-6489cee1c35d
function generate_queries(queries)
	function query_range()
		query_length = rand(query_length_range)
		query_start = rand(1:length(sequence) - query_length)
		query_end = query_start + query_length
		query_start:query_end
	end
	
	
	ranges = [query_range() for _ in 1:queries]
	samples = [LongDNA{2}(sample_range((@view sequence[range]), (@view  probabilites[range]))) for range in ranges]
	
	for sample in samples
		mutate!(sample)
	end
	
	(samples, ranges)
end

# ╔═╡ 50745be3-8539-48fa-9fd6-4d184ec30862
begin
	Random.seed!(0)
	probs = rand(0.8:0.01:1, word_size)
	queries, ranges = generate_queries(1)
	query = queries[1]
	range = ranges[1]
	subject = sequence[first(range)-50:last(range)+50]
end

# ╔═╡ 497ed7bc-d724-4d77-96be-036b6131675d
# ╠═╡ disabled = true
#=╠═╡
function seed_extension(query, seed_index, sequence_index)
	seed_start = seed_index
	seed_end = seed_index + word_size - 1
	q = query[seed_index:seed_index + word_size - 1]
	score = prob_word_score(q, subject, probs)
	
	while score > S_score_cutoff
		left_gap = query[seed_index:seed_index + word_size]
		pushfirst!(left_gap, DNA_Gap)
		
		right_gap = query[seed_index - 1:seed_index + word_size - 1]
		push!(right_gap, DNA_Gap)
		
		both_gap = query[seed_index:seed_index + word_size - 1]
		push!(both_gap, DNA_Gap)
		pushfirst!(both_gap, DNA_Gap)
		
		if seed_start != 1
			seed_start -= 1
		end
		if seed_end != length(query)
			seed_end += 1
		end
		
		no_gap = query[seed_index:seed_index + word_size - 1]
	end
end
  ╠═╡ =#

# ╔═╡ 91b9e822-1982-4d80-932c-c66ac2cb8a06
#=╠═╡
@time PBLAST_max3(query)
  ╠═╡ =#

# ╔═╡ 8d83a59b-962f-4546-957e-b1e06ceb3c41
#=╠═╡
@time PBLAST_max2(query)
  ╠═╡ =#

# ╔═╡ 9b7d489f-959c-4b9c-ba81-1720abf2b85f
function accuracy_test(seed, query_n)
	Random.seed!(seed)
	queries, ranges = generate_queries(query_n)
	matches_ = PBLAST_list(queries, results_n = 10)
	max_accuracy = proportionmap(matches(ranges, (first.(matches_))))
	#top_10_accuracy = proportionmap(map((r, matches)->r in [m.range for m in matches], ranges, matches))
	#(max_accuracy, top_10_accuracy)
end

# ╔═╡ 24aba1f1-186f-4a11-ad23-5dd6477e881d
accuracy_test(0, 100)

# ╔═╡ fd65ecce-d273-4551-8f77-213ea12d386a
begin
	Random.seed!(0)
	generate_queries(1)
end

# ╔═╡ af40cee4-e0eb-4764-b658-a2ff4d71d95e
# ╠═╡ disabled = true
#=╠═╡
@benchmark generate_queries(100) samples=50
  ╠═╡ =#

# ╔═╡ 28dbfc4b-63ea-4b10-822b-6e9e35dcfccc
# ╠═╡ disabled = true
#=╠═╡
@time PBLAST_list(queries, results_n = 1)
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
BioAlignments = "00701ae9-d1dc-5365-b64a-a3a3ebf5695e"
BioSequences = "7e6ae17a-c86d-528c-b3b9-7f778a29fe59"
Chain = "8be319e6-bccf-4806-a6f7-6fae938471bc"
Dictionaries = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
FASTX = "c2308a5c-f048-11e8-3e8a-31650f418d12"
Kmers = "445028e4-d31f-4f27-89ad-17affd83fc22"
Lazy = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
ProfileSVG = "132c30aa-f267-4189-9183-c8a63c7e05e6"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
Unzip = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"

[compat]
BenchmarkTools = "~1.3.2"
BioAlignments = "~3.0.0"
BioSequences = "~3.1.2"
Chain = "~0.5.0"
Dictionaries = "~0.3.25"
FASTX = "~2.0.0"
Kmers = "~0.1.0"
Lazy = "~0.15.1"
Parameters = "~0.12.3"
ProfileSVG = "~0.2.1"
StatsBase = "~0.33.21"
Unzip = "~0.2.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "2c06c76b7bc3e5031bbf19ce4ac157e7250981ab"

[[deps.AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["Printf", "ScanByte", "TranscodingStreams"]
git-tree-sha1 = "d50976f217489ce799e366d9561d56a98a30d7fe"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "0.8.2"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "d9a9701b899b30332bbcb3e1679c41cce81fb0e8"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.2"

[[deps.BioAlignments]]
deps = ["BioGenerics", "BioSequences", "BioSymbols", "IntervalTrees", "LinearAlgebra"]
git-tree-sha1 = "707e7e02cfb91a2b09e8f89df0204603b026df88"
uuid = "00701ae9-d1dc-5365-b64a-a3a3ebf5695e"
version = "3.0.0"

[[deps.BioGenerics]]
deps = ["TranscodingStreams"]
git-tree-sha1 = "0b581906418b93231d391b5dd78831fdc2da0c82"
uuid = "47718e42-2ac5-11e9-14af-e5595289c2ea"
version = "0.1.2"

[[deps.BioSequences]]
deps = ["BioSymbols", "Random", "SnoopPrecompile", "Twiddle"]
git-tree-sha1 = "e67b5446b44d96595ededa2b32ab0f2f6b3d2997"
uuid = "7e6ae17a-c86d-528c-b3b9-7f778a29fe59"
version = "3.1.2"

[[deps.BioSymbols]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "2052c3ec7c41b69efa0e9ff7e2734aa6658d4c40"
uuid = "3c28c6f8-a34d-59c4-9654-267d177fcfa9"
version = "5.1.2"

[[deps.Chain]]
git-tree-sha1 = "8c4920235f6c561e401dfe569beb8b924adad003"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.5.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e7ff6cadf743c098e08fca25c91103ee4303c9bb"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "00a2cccc7f098ff3b66806862d275ca3db9e6e5a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.5.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.DataAPI]]
git-tree-sha1 = "e08915633fcb3ea83bf9d6126292e5bc5c739922"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.13.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Dictionaries]]
deps = ["Indexing", "Random", "Serialization"]
git-tree-sha1 = "e82c3c97b5b4ec111f3c1b55228cebc7510525a2"
uuid = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
version = "0.3.25"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "c36550cb29cbe373e95b3f40486b9a4148f89ffd"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.2"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FASTX]]
deps = ["Automa", "BioGenerics", "BioSequences", "ScanByte", "StringViews", "TranscodingStreams"]
git-tree-sha1 = "9c72011edd523a83bf00276c4697a3019cea2257"
uuid = "c2308a5c-f048-11e8-3e8a-31650f418d12"
version = "2.0.0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.FlameGraphs]]
deps = ["AbstractTrees", "Colors", "FileIO", "FixedPointNumbers", "IndirectArrays", "LeftChildRightSiblingTrees", "Profile"]
git-tree-sha1 = "d9eee53657f6a13ee51120337f98684c9c702264"
uuid = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"
version = "0.2.10"

[[deps.Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalTrees]]
deps = ["InteractiveUtils", "Profile", "Random", "Test"]
git-tree-sha1 = "6c9fcd87677231ae293f6806fad928c216ab6658"
uuid = "524e6230-43b7-53ae-be76-1e9e4d08d11b"
version = "1.0.0"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.Kmers]]
deps = ["BioSequences"]
git-tree-sha1 = "d9d7b903c44f348f113659bf89e22d72a26fa048"
uuid = "445028e4-d31f-4f27-89ad-17affd83fc22"
version = "0.1.0"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LeftChildRightSiblingTrees]]
deps = ["AbstractTrees"]
git-tree-sha1 = "b864cb409e8e445688bc478ef87c0afe4f6d1f8d"
uuid = "1d6d02ad-be62-4b6b-8a6d-2f90e265016e"
version = "0.1.3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "946607f84feb96220f480e0422d3484c49c00239"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.19"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "b64719e8b4504983c7fca6cc9db3ebc8acc2a4d6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.ProfileSVG]]
deps = ["Colors", "FlameGraphs", "Profile", "UUIDs"]
git-tree-sha1 = "e4df82a5dadc26736f106f8d7fc97c42cc6c91ae"
uuid = "132c30aa-f267-4189-9183-c8a63c7e05e6"
version = "0.2.1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
git-tree-sha1 = "bc12e315740f3a36a6db85fa2c0212a848bd239e"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.2"

[[deps.ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "2436b15f376005e8790e318329560dcc67188e84"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.3.3"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StringViews]]
git-tree-sha1 = "609585ed628a4cd46f4c142762be37f5ced5dc7d"
uuid = "354b36f9-a18e-4713-926e-db85100087ba"
version = "1.0.3"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "e4bdc63f5c6d62e80eb1c0043fcc0360d5950ff7"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.10"

[[deps.Twiddle]]
git-tree-sha1 = "29509c4862bfb5da9e76eb6937125ab93986270a"
uuid = "7200193e-83a8-5a55-b20d-5d36d44a0795"
version = "1.1.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═f1d327bf-42ce-480e-b98c-42ac5ba67235
# ╠═4a7dc282-7355-4427-ad0a-efc1c3f9faad
# ╠═74f5561d-2039-40f0-88a2-b19fe555c3a2
# ╟─9e7f6e3a-b646-4fe4-bf8b-26159f816505
# ╠═4af62ee8-f649-400f-8182-1ac42ded6ae5
# ╠═4d4578e4-4c67-45bd-a360-70193aae5d5e
# ╠═5625822c-2030-4303-88d1-0c5e95cae278
# ╟─4c4979f3-dd90-4744-8c55-8eca1f9ccbea
# ╟─6db80ce9-4c47-4672-b0c6-7b81424ceb71
# ╟─57c5982a-1df9-4cc2-9d46-537cd31e2d68
# ╠═f4857a38-02bb-4c02-89d4-aa4906229b3e
# ╠═9edfc292-0aa5-4e07-8381-d7d02ec1a09d
# ╠═62d5b131-fd26-4225-8420-cd8b42a7200b
# ╠═ccc968fa-f75f-4a30-90d1-cf0c97288ae0
# ╠═c8948347-6cd9-41ab-aaf3-2e108ce02967
# ╟─eb97c760-a32f-42d5-8c7b-f9138b576fe2
# ╟─61b0eb60-f1c4-4067-9d70-ee181dd64a36
# ╠═d43926d5-d305-4e12-bfc3-d72d6aa8c000
# ╠═c0158950-9a2e-42d3-8d82-3249cad660a2
# ╟─3220651e-c49b-4cc6-b76d-26790451affc
# ╠═c485200c-8637-45f9-85e6-df399758a09e
# ╟─b5d3d606-91d6-454c-846b-81c744914042
# ╟─699cfd94-5271-4aaa-8bae-5772324f58a3
# ╟─3559f2e9-9a00-44f7-bc7e-1c05eab7af27
# ╟─ef8a4d03-e6fb-4e8e-9a8c-54f0ac2d8c7f
# ╟─fa9cfa41-26d6-4817-a063-ca4e46b91f44
# ╟─fdff2519-1b53-45a4-ad00-8bb4a382bf9d
# ╠═8b7aab0c-3690-4549-b840-e704ad1b70e4
# ╠═a1dc1cb0-58c5-44cf-8256-33d19c97b231
# ╠═58be55eb-d107-4d20-b2f2-a084f83ae338
# ╠═554206ee-5a38-47c3-9c68-d92fafbd6401
# ╟─0100a578-9f44-41a0-b0f0-dc2332d42325
# ╟─373d13c9-c11f-4a03-8ea5-5aef62471c4c
# ╟─497ed7bc-d724-4d77-96be-036b6131675d
# ╟─99b513ff-9928-4aba-8818-502c9136cbb0
# ╟─22bca925-c356-4a4c-a963-d2ffb42af719
# ╟─7e3d62bb-0892-4e5a-991c-0fbbfd7b2ae9
# ╟─ae3725c3-0f0c-431c-8217-ca55c66aac05
# ╠═50745be3-8539-48fa-9fd6-4d184ec30862
# ╠═91b9e822-1982-4d80-932c-c66ac2cb8a06
# ╠═8d83a59b-962f-4546-957e-b1e06ceb3c41
# ╠═9ecd59ea-174b-4ad0-b943-d05566088bc3
# ╠═f7470de9-4a5a-4395-b2c5-ccf91734824b
# ╠═623e3c50-8d86-418f-98c2-7a9c9da769bd
# ╠═0035cbb3-6f3b-4e04-aa2c-d125b8baffbe
# ╠═5153c32d-2881-4461-9077-3476ced37798
# ╠═69559e82-65ee-4986-95f7-eead25f83e56
# ╟─ea42ca53-d31f-4370-a350-1f1f9fde8f41
# ╟─6e2ed658-bc97-4749-8b84-243c534c7d34
# ╠═4870eeeb-cd46-41d4-98e7-f87af65b3f21
# ╟─59a677a6-c2f4-4921-bcc8-dbe27242fc44
# ╟─ed31ab1c-70ea-4f6d-9824-f4cb6f8dc48d
# ╟─d13c8a85-e7d5-4a81-a70a-57b8ab6206a0
# ╠═88fe623f-6163-4e8e-8e34-5ce967c84ae2
# ╠═8fa9fd5c-ebd9-46c5-9695-471aef7b156c
# ╠═0c6c200f-732d-4831-95f4-5f977e7a2b48
# ╠═81f8a645-53ff-4e67-b6cf-25e3343c4b04
# ╠═e92a87b0-d37a-4b1b-9740-fe89a963cb16
# ╠═4980de9a-b1c8-49b0-9db2-bc38f76b13ee
# ╠═cf34efa1-87ab-4094-a667-cf679f6edecf
# ╠═9b7d489f-959c-4b9c-ba81-1720abf2b85f
# ╠═ac9c68b3-0a7a-468d-9b84-9c4a668c3c32
# ╠═24aba1f1-186f-4a11-ad23-5dd6477e881d
# ╠═8006ffbb-e654-4aff-9666-6489cee1c35d
# ╠═5baedc9a-654b-4b96-8150-324ecf94c4ee
# ╠═fd65ecce-d273-4551-8f77-213ea12d386a
# ╠═af40cee4-e0eb-4764-b658-a2ff4d71d95e
# ╟─28dbfc4b-63ea-4b10-822b-6e9e35dcfccc
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
