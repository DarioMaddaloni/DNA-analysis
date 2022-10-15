""" Task 1 """
### Dario Maddaloni
def kmer_hist(S,k): #S is a sequence and n is its length. k < n is the length of the strings we are interested in
    ############ Data collection ##############
    kmer=S[0:k] #kmer will be the sequece of characters we are interested in (ha s length k)

    my_dict={kmer:1} # my_dict will be a dictionary where the key is the kmer and the item is the frequency
    max_freq=1 # we set max_freq to 1 because we have already initialized my_dict with the first k-mer. Hence, we can state that kmer has frequency 1 until now
    
    for char in S[k:]: # we start from k position because the first kmer take up k characters (but the indexes of the sequence start from 0)
        kmer=kmer[1:]+char # we update the kmer to the following one
        if kmer in my_dict.keys(): # it costs a lot: it is a search in a sequence
            my_dict[kmer]+=1
            if max_freq < my_dict[kmer]:
                max_freq=my_dict[kmer]
        else: # if there are few repetition the "cycle for" costs less
            my_dict[kmer]=1
    
    
    ############## Data organization ##############
    # Now we start to reorganize the dectionary for the output
    h=[0 for i in range(max_freq+1)] # max_freq+1 beacuse we want also that h[0]=0
    my_second_dict={0:{""}} #keys are the frequencies, items are kmers. It is trivial to see that there are exacly 0 kmer (in the set of all the kmer that appear in the sequence S) that repeat themselves 0 times
    for kmer,freq in my_dict.items(): # freq represents the frequency of the item kmer in my first dictionary
        h[freq]+=1 # there is a kmer more with frequency freq
        if h[freq]==1: # if the statement returns true, this means that we are in the case that we have never found a kmer of frequency freq
            my_second_dict[freq]={kmer}
        else:
            my_second_dict[freq].add(kmer)
            
    
    mfkmers=list(my_second_dict[len(h)-1])
    
    return h, mfkmers;
# at the position i of the sequence h we have the number of substrings of S that have frequency i
# we have to note that h[0]=0 always and so the mfkmer is in the position len(h)-1 and the frequence is len(n)-1
# mfkmers contain the most frequent kmers. They can be more than one, so mfkmers is a sequence.



""" Task 3 """
def kmer_search(S,L):
    pos, freq = None, 0 # if I can't find any of the kmers in L, these are the values I will return 
    
    seq_pos=[None for i in range(len(L))]#is the sequence of all the positions
    seq_freq=[0 for i in range(len(L))]#is the sequence of all the frequencies
    
    k=len(L[0]) # k is the length of the k_mers. I assume that they have the same length
    string=S[0:k] # string is now the first kmer in the sequence S
    
    ############## Subsection1 ##############
    for j in range(len(S)-k-1): # j+k is the position of the character we want to append to the kmer we want to control in L. We use a cycle on len(S)-k instead of S[k:] because we will use j as the position we want to save in our seq_pos
        string=string[1:]+S[j+k] #I'll control the following kmer
        for i in range(len(L)): #this is my i-th check if the string S[j:j+k] is in L
            if string==L[i]:
                seq_freq[i]+=1
                if seq_pos[i]==None:
                    seq_pos[i]=j+1 #I have to take the position (+1, remember the index starts from 0) of the first character of the string
                break # we can save some rounds
    #Now, I have seq_freq with all the frequencies of L and seq_pos with the first position that L[i] occurs
    
    ############## maximum_pos and maximum_freq ##############
    #Now I want to find maximum_pos and maximum_freq
    for i in range(len(seq_pos)): # I don't want to inizialize maximum_pos with None so I do this "for cycle"
        if seq_pos[i]!=None:
            maximum_freq=seq_freq[i] # there could be some None
            maximum_pos=seq_pos[i] # it will be the position of the most frequent kmer in the list L and the first that appear in S
            break
            
    for i in range(len(seq_freq)):
        if seq_freq[i]>=maximum_freq and seq_pos[i]!=None: #if the frequence of L[i] is higher or equal than the maximum of the frequencies of L[:i]. Moreover the second statement means that I actually found that kmer
            maximum_pos=seq_pos[i]
            maximum_freq=seq_freq[i] #we are not interested in which of the two relations seq_freq and maximum_freq are: indeed, we are interested only in the fact that this assignment is non-decreasing
            
    return maximum_pos, maximum_freq
