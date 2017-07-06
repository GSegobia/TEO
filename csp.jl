function hamming(sentence1, sentence2)
    distance = 0
    for i in 1:length(sentence1)
        if sentence1[i] != sentence2[i]
            distance += 1
        end
    end
    return distance
end

function build(n, string_length, data_collection, α)
    solution = Array{Char}(string_length)
    candidate_list = Array{Int64}(n)
    distance_list = Array{Int64}(n)
    rcl = []
    iMax = 0
    distance_matrix = zeros(length(data_collection),length(data_collection))
    for i in 1:n
        candidate_list[i] = 0
        distance_list[i] = 0
    end
    for i in 1:string_length
        max_distance = 0
        min_distance = 999999
        for j in 1:n
            for k in 1:n
                # println(data_collection[i][j])
                # println(data_collection[i][k])
                # distance_matrix[1,1]=0
                if data_collection[j][i] != data_collection[k][i]
                    distance_matrix[k,j] = 1
                else
                    distance_matrix[k,j] = 0
                end
                if (distance_matrix[j,k] + distance_list[k]) > candidate_list[k]
                    candidate_list[j] = distance_matrix[j,k] + distance_list[k]
                end
            end
            if candidate_list[j] > max_distance
                max_distance = candidate_list[j]
            end
            if candidate_list[j] < min_distance
                min_distance = candidate_list[j]
            end
        end
        for j in 1:n
            if min_distance + ((max_distance - min_distance) * α) >= candidate_list[j]
                push!(rcl, data_collection[j][i])
            end
        end
        # println(candidate_list)
        # println(rcl)
        solution[i] = rcl[rand(1:length(rcl))]
        iMax = -1
        for j in 1:n
            #println(hamming(solution,String(data_collection[j])))
            distance_list[j] = hamming(solution,data_collection[j])
            if distance_list[j] >= iMax
                iMax = distance_list[j]
            end
            candidate_list[j] = 0
        end
        empty!(rcl)
    end
    return iMax,solution,distance_list
end

function local_search(iMax, solution, distance_list, data_collection, n, string_length, alphabet)
    solution_l = copy(solution)
    for i in 1:string_length
        iMax_l = -1
        alphabet_l = copy(alphabet)
        deleteat!(alphabet_l, find(x -> x == solution[i], alphabet_l))
        random_value = alphabet_l[rand(1:length(alphabet_l))]
        solution_l[i] = random_value
        for j in 1:n
            distance_list[j] = hamming(solution_l, data_collection[j])
            if distance_list[j] >= iMax_l
                iMax_l = distance_list[j]
            end
        end
        if iMax_l < iMax
            solution[i] = solution_l[i]
            iMax = iMax_l
        else
            solution_l[i] = solution[i]
        end
    end

    return iMax, solution
end

alphabet = ['A', 'C', 'T', 'G']

dna_basis = []

# dna_basis = ["CCAGCTGCATCACAGGAGGCCAGCGAGCAGGTCTGTTCCAAGGGCCTTCGAGCCAGTCTG",
#                   "AGACCCGCCGGGAGGCGGAGGACCTGCAGGGTGAGCCCCACCGCCCCTCCGTGCCCCCGC",
#                   "GAGGTGAAGGACGTCCTTCCCCAGGAGCCGGTGAGAAGCGCAGTCGGGGGCACGGGGATG",
#                   "GGGCTGCGTTGCTGGTCACATTCCTGGCAGGTATGGGGCGGGGCTTGCTCGGTTTTCCCC",
#                   "GCTCAGCCCCCAGGTCACCCAGGAACTGACGTGAGTGTCCCCATCCCGGCCCTTGACCCT",
#                   "CAGACTGGGTGGACAACAAAACCTTCAGCGGTAAGAGAGGGCCAAGCTCAGAGACCACAG"]

open("Data/ins.data") do file
    while !eof(file)
        #push!(dna_basis, strip(readline(file), '\n')) # Rodar essa linha no linux, pq o windows eh uma máquina de escrever
        push!(dna_basis, strip(readline(file), ['\n', '\r'])) #Eu odeio o windows, maldito Bill Gates
    end
end

N = length(dna_basis)
M = length(dna_basis[1])

println("Número de bases: $(N)")
println("Tamanho das bases: $(M)")

results = build(N, M, dna_basis, 0.8)
println("Menor distância máxima fase de construção: $(results[1])")
println("Sequência encontrada fase de contrução: $(String(results[2]))")

results_build = [copy(results[1]), String(results[2])]

new_results = local_search(results[1], results[2], results[3], dna_basis, N, M, alphabet)
println("Menor distância máxima busca local: $(new_results[1])")
println("Sequência encontrada busca local: $(String(new_results[2]))")
