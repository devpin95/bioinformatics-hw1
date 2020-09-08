gene = 'CXCR6'
fullstr = ''
codons = []

aa_dict = {
    'u': {
        'u': {'u': 'F', 'c': 'F', 'a': 'L', 'g': 'L'},
        'c': {'u': 'S', 'c': 'S', 'a': 'S', 'g': 'S'},
        'a': {'u': 'Y', 'c': 'Y', 'a': 'S', 'g': 'S'},
        'g': {'u': 'C', 'c': 'C', 'a': '*', 'g': 'W'},
    },
    'c': {
        'u': {'u': 'L', 'c': 'L', 'a': 'L', 'g': 'L'},
        'c': {'u': 'P', 'c': 'P', 'a': 'P', 'g': 'P'},
        'a': {'u': 'H', 'c': 'H', 'a': 'Q', 'g': 'Q'},
        'g': {'u': 'R', 'c': 'R', 'a': 'R', 'g': 'R'},
    },
    'a': {
        'u': {'u': 'I', 'c': 'I', 'a': 'I', 'g': '^M'},
        'c': {'u': 'T', 'c': 'T', 'a': 'T', 'g': 'T'},
        'a': {'u': 'N', 'c': 'N', 'a': 'K', 'g': 'K'},
        'g': {'u': 'S', 'c': 'S', 'a': 'R', 'g': 'R'},
    },
    'g': {
        'u': {'u': 'V', 'c': 'V', 'a': 'V', 'g': 'V'},
        'c': {'u': 'A', 'c': 'A', 'a': 'A', 'g': 'A'},
        'a': {'u': 'D', 'c': 'D', 'a': 'E', 'g': 'E'},
        'g': {'u': 'G', 'c': 'G', 'a': 'G', 'g': 'G'},
    }
}


def codon_to_aa(codonstr):
    return aa_dict[codonstr[0]][codonstr[1]][codonstr[2]]


with open('genes/' + gene) as gcoding:
    for x in gcoding:
        line = x.split(' ')
        line.pop(0)
        line[-1] = line[-1].rstrip('\r\n')
        for seq in line:
            fullstr = fullstr + seq

for i in range(0, len(fullstr), 3):
    codon = fullstr[i]

    if i+1 < len(fullstr):
        codon = codon + fullstr[i+1]

    if i+2 < len(fullstr):
        codon = codon + fullstr[i+2]

    codon = codon.replace('t', 'u')

    codons.append(codon)

start_codon = 0
for i in range(0, len(codons)):
    if codons[i] == 'aug':
        start_codon = i
        break

stop_codon = 0
for i in range(0, len(codons)):
    if codon_to_aa(codons[i]) == '*' and i > start_codon:
        stop_codon = i
        break

print('Start Codon index: ' + str(start_codon))
print('End Codon index: ' + str(stop_codon))

print('\nDNA (start codon to end codon)')
print(codons)
print(codons[start_codon] + '-' + codons[start_codon+1] + '-' + codons[start_codon+2] + '-' + codons[start_codon+3] + '-...-' + codons[stop_codon-3] + '-' + codons[stop_codon-2] + '-' + codons[stop_codon-1] + '-' + codons[stop_codon])
print('M-' + codon_to_aa(codons[start_codon+1]) + '-' + codon_to_aa(codons[start_codon+2]) + '-' + codon_to_aa(codons[start_codon+3]) + '-...-' + codon_to_aa(codons[stop_codon-3]) + '-' + codon_to_aa(codons[stop_codon-2]) + '-' + codon_to_aa(codons[stop_codon-1]) + '-*')

translation = ''
for codon in codons:
    if len(codon) == 3:
        translation = translation + codon_to_aa(codon)

first_start_codon = translation.find('^')
first_stop_codon = translation.find('*')

translation_clean = translation[first_start_codon : first_stop_codon].replace('^', '')


print('\nTranslation')
for i in range(0, len(translation_clean)):
    print(translation_clean[i], end='')
    if i % 50 == 0 and not i == 0:
        print('')


# MEKARPLWANSLQFVFACISYAVGLGNVWRFPYLCQMYGGGSFLVPYIIML
# IVEGMPLLYLELAVGQRMRQGSIGAWRTISPYLSGVGVASVVVSFFLSMY
# YNVINAWAFWYLFHSFQDPLPWSVCPLNGNHTGYDEECEKASSTQYFWYR
# KTLNISPSLQENGGVQWEPALCLLLAWLVVYLCILRGTESTGKVVYFTAS
# LPYCVLIIYLIRGLTLHGATNGLMYMFTPKIEQLANPKAWINAATQIFFS
# LGLGFGSLIAFASYNEPSNNCQKHAIIVSLINSFTSIFASIVTFSIYGFK
# ATFNYENCLKKVSLLLTNTFDLEDGFLTASNLEQVKGYLASAYPSKYSEM
# FPQIKNCSLESELDTAVQGTGLAFIVYTEAIKNMEVSQLWSVLYFFMLLM
# LGIGSMLGNTAAILTPLTDSKIISSHLPKEAISGLVCLVNCAIGMVFTME
# AGNYWFDIFNDYAATLSLLLIVLVETIAVCYVYGLRRFESDLKAMTGRAV
# SWYWKVMWAGVSPLLIVSLFVFYLSDYILTGTLKYQAWDASQGQLVTKDY
# PAYALAVIGLLVASSTMCIPLAALGTFVQRRLKRGDADPVA