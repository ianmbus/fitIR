


import cPickle as pickle

# class Fruits: pass
# 
# banana = Fruits()
# 
# banana.color = 'yellow'
# banana.value = 30
# 
# 
# pickle.dump(banana, open('test.pck','w'), pickle.HIGHEST_PROTOCOL)
# 


athing = pickle.load(open('test.pck','r'))

print athing
print athing.color
