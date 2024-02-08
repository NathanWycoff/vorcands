#from bnn-bo import distillation
#from bnn_bo import distillation
#from bnn_bo.test_functions.distillation import KnowledgeDistillation

import numpy as np
import os
os.environ["SDL_VIDEODRIVER"] = "dummy"

#exec(open('../bnn-bo/test_functions/distillation.py').read())
exec(open('./python/meta_module.py').read())

exec(open('./python/lunar_lander.py').read())
prob_lunar = LunarLanderProblem()
lunar = lambda x: float(-prob_lunar.get_reward(x))
#m = 12
#x = np.random.uniform(size=m)
#y = lunar(x)

exec(open('./python/pdes.py').read())
prob_pde = PDEVar()
#xc = (prob_pde.bounds[1,:]-prob_pde.bounds[0,:]) * x + prob_pde.bounds[0,:]
Ppde = lambda x: float(prob_pde.evaluate_true([(prob_pde.bounds[1,:]-prob_pde.bounds[0,:]) * torch.tensor(x) + prob_pde.bounds[0,:]])[0])
#Ppde = lambda x: prob_pde.evaluate_true([x])[0]
#m = 4
#x = np.random.uniform(size=m)
#Ppde(x)

