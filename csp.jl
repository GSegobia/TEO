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
    return iMax,solution,length(solution)
end


dna_basis = []

dna_basis_test = ["CCAGCTGCATCACAGGAGGCCAGCGAGCAGGTCTGTTCCAAGGGCCTTCGAGCCAGTCTG",
                  "AGACCCGCCGGGAGGCGGAGGACCTGCAGGGTGAGCCCCACCGCCCCTCCGTGCCCCCGC",
                  "GAGGTGAAGGACGTCCTTCCCCAGGAGCCGGTGAGAAGCGCAGTCGGGGGCACGGGGATG",
                  "GGGCTGCGTTGCTGGTCACATTCCTGGCAGGTATGGGGCGGGGCTTGCTCGGTTTTCCCC",
                  "GCTCAGCCCCCAGGTCACCCAGGAACTGACGTGAGTGTCCCCATCCCGGCCCTTGACCCT",
                  "CAGACTGGGTGGACAACAAAACCTTCAGCGGTAAGAGAGGGCCAAGCTCAGAGACCACAG"]

open("Data/splice.data") do file
    while !eof(file)
        push!(dna_basis, strip(readline(file), '\n'))
    end
end

N = length(dna_basis_test)
M = length(dna_basis_test[1])

println(N)
println(M)

results = build(N, M, dna_basis_test, 0.9)

println(results[1])
println(String(results[2]))
println(results[3])
