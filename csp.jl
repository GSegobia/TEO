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
        solution[i] = rcl[rand(1:length(rcl))]
        iMax = -1
        for j in 1:n
            distance_list[j] = hamming(solution,data_collection[j])
            if distance_list[j] >= iMax
                iMax = distance_list[j]
            end
            candidate_list[j] = 0
        end
        empty!(rcl)
    end
    return iMax, solution
end

#Determinístico e Best Improvement
function local_search_grasp2(iMax, solution, data_collection, number_of_strings, string_length, alphabet)
    distance_list = Array{Int64}(number_of_strings)
    stop_condition = true
    while(stop_condition)
        random_value = alphabet[rand(1:length(alphabet))]
        solution_l = copy(solution)
        iMax_l = iMax
        for i in 1 : string_length
            neighboor = copy(solution)
            iMax_neighboor = -1
            neighboor[i] = random_value
            for j in 1 : number_of_strings
                distance_list[j] = hamming(neighboor, data_collection[j])
                if distance_list[j] > iMax_neighboor
                    iMax_neighboor = distance_list[j]
                end
            end
            if iMax_l > iMax_neighboor
                iMax_l = iMax_neighboor
                solution_l[i] = neighboor[i]
            end
        end
        if iMax_l < iMax
            iMax = iMax_l
            solution = copy(solution_l)
        else
            stop_condition = false
        end
    end
    return iMax, solution
end

#Aleatório e First Improvement
function local_search_grasp3(iMax, solution, data_collection, number_of_strings, string_length, alphabet)
    distance_list = Array{Int64}(number_of_strings)
    stop_condition = true
    count = 0
    while(stop_condition)
        random_alphabet_value = alphabet[rand(1:length(alphabet))]
        iMax_control = iMax
        for i in 1 : string_length
            random_position_value = rand(1:string_length)
            solution_l = copy(solution)
            solution_l[random_position_value] = random_alphabet_value
            iMax_l = -1
            for j in 1 : number_of_strings
                distance_list[j] = hamming(solution_l, data_collection[j])
                if distance_list[j] > iMax_l
                    iMax_l = distance_list[j]
                end
            end
            if iMax_l < iMax
                solution[random_position_value] = solution_l[random_position_value]
                iMax = iMax_l
                break
            end
        end
        if iMax == iMax_control
            count = count + 1
            if count == 10
                stop_condition = false
            end
        end
    end
    return iMax, solution
end

function local_search(iMax, solution, data_collection, number_of_strings, string_length, alphabet)
    distance_list = Array{Int64}(number_of_strings)
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

open("instancias/n10m1000tai4.ins") do file
    readline(file)
    readline(file)
    readline(file)
    readline(file)
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
#println("Sequência encontrada fase de contrução: $(String(results[2]))")

results_build = [copy(results[1]), String(results[2])]

new_results = local_search_grasp3(results[1], results[2], dna_basis, N, M, alphabet)
println("Menor distância máxima busca local: $(new_results[1])")
#println("Sequência encontrada busca local: $(String(new_results[2]))")
