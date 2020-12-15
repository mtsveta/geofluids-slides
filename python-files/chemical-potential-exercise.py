mu_case1 = [-239.614, -390.347, -610.019, -19.942]
mu_case2 = [-239.614, -390.347, -717.639, -41.454]
mu_case3 = [-239.614, -389.185, -606.776, -16.699]

def reaction_direction(mu):
    diff = (mu[0] + mu[1]) - (mu[2] + mu[3])
    if diff == 0.0:
        return 'equilibrium'
    elif diff > 0.0:
        return 'right'
    else:
        return 'left'

print( 'Case I:  ', reaction_direction(mu_case1) )
print( 'Case II: ', reaction_direction(mu_case2) )
print( 'Case III:', reaction_direction(mu_case3) )