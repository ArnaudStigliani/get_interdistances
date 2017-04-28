import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker
import numpy as np
import sys
from matplotlib import gridspec

def divide(a, b):
    if b == 0:
        return np.nan
    else: 
        return a/b

def negative_sets_histo (Interdistance_maxValue,relative_DR,relative_DR_neg,relative_ER,relative_ER_neg,relative_IR,relative_IR_neg,factorTranscription,threshold,rate,command) :
	fig = plt.figure(1,figsize= (18,10))
	indexes1 = range(Interdistance_maxValue + 1)
	width = 1
	ax1 = fig.add_subplot(1,4,1)
	ax1.set_xlabel("", size = 16)
	indexes1 = np.arange(3)
	plt.xticks(indexes1 + width * 0.5, ('DR', 'ER', 'IR'))
	ax1.axis([0, 3, 0, 2])
	ax1.bar(indexes1, rate , width , color = 'green')
	ax1.set_ylabel('$\sum DR$n_+, ER$n_+, IR$n_+ / $\sum DR$n_-, ER$n_-, IR$n_- , respectively', color = 'green', size = 16)
	print("relative_DR : ",relative_DR)
	print("relative_DR_neg : ",relative_DR_neg)
	for a,b in zip(relative_DR,relative_DR_neg) :
		ax1 = fig.add_subplot(1,4,2)
		ax1.set_xlabel("base pairs between "+factorTranscription+" direct repeats", size = 16)
		indexes1 = range(Interdistance_maxValue + 1)
		ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		print("len(a) : ",len(a))
		print("len(b) : ",len(b))
		ax1.bar(indexes1, map(divide, a, b) , width , color = 'cornflowerblue')
		ax1.set_ylabel(' ', color = 'cornflowerblue', size = 16)
		ax2 = ax1.twinx()
		ax2.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		ax2.bar(indexes1,[x * float(10) for x in a] , width , color = 'b')
		ax2.bar(indexes1,[x * float(10) for x in b] , width, color = 'brown', alpha = 0.85)
		ax2.set_ylabel('DR$n_+$ frequence / D$Rn_-$ frequence', color = 'cornflowerblue', size = 16)
		ybox1 = TextArea("DR$n_+$ frequence (X10)", textprops = dict(color = "b", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
		ybox3 = TextArea("DR$n_-$ frequence (X10), ", textprops = dict(color = "brown", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
		ybox = VPacker(children = [ybox1, ybox3],align = "bottom", pad = 0, sep = 8)
		anchored_ybox = AnchoredOffsetbox(loc = 8, child = ybox, pad = 0., frameon = False, bbox_to_anchor = (-0.08, 0.3), 
				bbox_transform = ax2.transAxes, borderpad = 0.)
		ax2.add_artist(anchored_ybox)
		indexes1 = np.arange(Interdistance_maxValue + 1)
		print("indexes1 : ",indexes1)
		print("width : ",width)
		plt.xticks(indexes1 + width * 0.5 , indexes1)
		plt.text(2, 5, "D$Rn_+$ = DRn number in the bound set")
		plt.text(2, 4, "D$Rn_-$ = DRn number in the unbound set")
		plt.text(2, 3, "DRn frequence = DRn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$")
		plt.text(2, 6, sys.argv)
		l = plt.axhline(y = 1)
		ax1.legend()
		ax2.legend()
		
	for a, b,c in zip(relative_ER,relative_ER_neg,threshold)  :	
		ax1 = fig.add_subplot(1,4,3)
		ax1.set_xlabel("base pairs between "+factorTranscription+" everted repeats", size = 16)
		indexes1 = range(Interdistance_maxValue + 1)
		ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		ax1.bar(indexes1, map(divide, a, b) , width , color = 'cornflowerblue')
		ax1.set_ylabel(' ', color = 'cornflowerblue', size = 16)
		ax2 = ax1.twinx()
		ax2.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		ax2.bar(indexes1,[x * float(10) for x in a] , width , color = 'b')
		ax2.bar(indexes1,[x * float(10) for x in b] , width, color = 'brown', alpha = 0.85)
		ax2.set_ylabel('E$Rn_+$ frequence / D$Rn_-$ frequence', color = 'cornflowerblue', size = 16)
		ybox1 = TextArea("E$Rn_+$ frequence (X10)", textprops = dict(color = "b", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
		ybox3 = TextArea("E$Rn_-$ frequence (X10), ", textprops = dict(color = "brown", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
		ybox = VPacker(children = [ybox1, ybox3],align = "bottom", pad = 0, sep = 8)
		anchored_ybox = AnchoredOffsetbox(loc = 8, child = ybox, pad = 0., frameon = False, bbox_to_anchor = (-0.08, 0.3), 
				bbox_transform = ax2.transAxes, borderpad = 0.)
		ax2.add_artist(anchored_ybox)
		indexes1 = np.arange(Interdistance_maxValue + 1)
		plt.xticks(indexes1 + width * 0.5 , indexes1)
		plt.text(2, 5, "ER$n_+$ = ERn number in the bound set")
		plt.text(2, 4, "ER$n_-$ = ERn number in the unbound set")
		plt.text(2, 3, "ERn frequence = ERn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$")
		l = plt.axhline(y = 1)
		plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
		ax1.legend()
		ax2.legend()
	for a,b in zip(relative_IR,relative_IR_neg) :	
		ax1 = fig.add_subplot(1,4,4)
		ax1.set_xlabel("base pairs between "+factorTranscription+" inverted repeats", size = 16)
		indexes1 = range(Interdistance_maxValue + 1)
		ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		ax1.bar(indexes1, map(divide, a, b) , width , color = 'cornflowerblue')
		ax1.set_ylabel(' ', color = 'cornflowerblue', size = 16)
		ax2 = ax1.twinx()
		ax2.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		ax2.bar(indexes1,[x * float(10) for x in a] , width , color = 'b')
		ax2.bar(indexes1,[x * float(10) for x in b] , width, color = 'brown', alpha = 0.85)
		ax2.set_ylabel('I$Rn_+$ frequence / D$Rn_-$ frequence', color = 'cornflowerblue', size = 16)
		ybox1 = TextArea("I$Rn_+$ frequence (X10)", textprops = dict(color = "b", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
		ybox3 = TextArea("I$Rn_-$ frequence (X10), ", textprops = dict(color = "brown", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
		ybox = VPacker(children = [ybox1, ybox3],align = "bottom", pad = 0, sep = 8)
		anchored_ybox = AnchoredOffsetbox(loc = 8, child = ybox, pad = 0., frameon = False, bbox_to_anchor = (-0.08, 0.3), 
				bbox_transform = ax2.transAxes, borderpad = 0.)
		ax2.add_artist(anchored_ybox)
		indexes1 = np.arange(Interdistance_maxValue + 1)
		plt.xticks(indexes1 + width * 0.5 , indexes1)
		plt.text(2, 5, "IR$n_+$ = DRn number in the bound set")
		plt.text(2, 4, "IR$n_-$ = DRn number in the unbound set")
		plt.text(2, 3, "IRn frequence = IRn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$")
		l = plt.axhline(y = 1)
		ax1.legend()
		ax2.legend()
	
	plt.show()
		
def negative_sets_curve(Interdistance_maxValue,relative_DR,relative_DR_neg,relative_ER,relative_ER_neg,relative_IR,relative_IR_neg,factorTranscription,threshold,rate,command) :
	fig = plt.figure(1,figsize= (18,10))
	indexes1 = range(Interdistance_maxValue + 1)
	width = 1
	ax1 = fig.add_subplot(1,4,1)
	ax1.set_xlabel("", size = 16)
	indexes1 = np.arange(3)
	plt.xticks(indexes1 + width * 0.5, ('DR', 'ER', 'IR'))
	ax1.axis([0, 3, 0, 2])
	ax1.bar(indexes1, rate , width , color = 'green')
	ax1.set_ylabel('$\sum DR$n_+, ER$n_+, IR$n_+ / $\sum DR$n_-, ER$n_-, IR$n_- , respectively', color = 'green', size = 16)
	plt.text(2, 3, command)
	for a,b in zip(relative_DR,relative_DR_neg) :
		indexes1 = range(Interdistance_maxValue + 1)
		ax1 = fig.add_subplot(1,4,2)
		ax1.set_xlabel("base pairs between "+factorTranscription+" direct repeats", size = 16)
		ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		plt.plot(indexes1, map(divide, a, b), lw=2,label=r"")
		ax1.set_ylabel("DR$n_+ frequence / DR$n_- frequence", size = 16)
		l = plt.axhline(y = 1)
		plt.legend()
	for a, b,c in zip(relative_ER,relative_ER_neg,threshold)  :	
		indexes1 = range(Interdistance_maxValue + 1)
		ax1 = fig.add_subplot(1,4,3)
		ax1.set_xlabel("base pairs between "+factorTranscription+" everted repeats", size = 16)
		ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		plt.plot(indexes1, map(divide, a, b), lw=2, label= r"threshold : "+str(c))
		ax1.set_ylabel("$ERn_+$ frequence / $ERn_-$ frequence", size = 16)
		l = plt.axhline(y = 1)
		plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
		plt.legend()
	for a,b in zip(relative_IR,relative_IR_neg) :
		indexes1 = range(Interdistance_maxValue + 1)
		ax1 = fig.add_subplot(1,4,4)
		ax1.set_xlabel("base pairs between "+factorTranscription+" inverted repeats", size = 16)
		ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		plt.plot(indexes1, map(divide, a, b), lw=2,label=r"")
		ax1.set_ylabel("$IRn_+$ frequence / $IRn_-$ frequence", size = 16)
		l = plt.axhline(y = 1)
		plt.legend()
	plt.legend()
	plt.show()
	
def negative_sets_points(Interdistance_maxValue,relative_DR,relative_DR_neg,relative_ER,relative_ER_neg,relative_IR,relative_IR_neg,factorTranscription,threshold,rate,command) :	
	fig = plt.figure(1,figsize= (18,10))
	gs = gridspec.GridSpec(1, 4, width_ratios=[1, 2,2,2]) 
	indexes1 = range(Interdistance_maxValue + 1)
	width = 1
	ax0 = plt.subplot(gs[0])
	ax0.set_xlabel("", size = 16)
	indexes1 = np.arange(3)
	plt.xticks(indexes1 + width * 0.5, ('DR', 'ER', 'IR'))
	ax0.axis([0, 3, 0, 2.5])
	print("rate : ",rate)
	ax0.bar(indexes1, rate , width , color = 'w')
	ax0.set_ylabel(r'$\sum$' + ' DRn+, ERn+, IRn+ / ' + r'$\sum$' + 'DRn-, ERn-, IRn- , respectively with a threshold of -12', size = 16)
	points = ['bo','go','ro','co','mo','yo','ko']
	turn = 0
	plt.text(-3.6, -0.2, command)
	
	maxi = []
	for a,b in zip(relative_DR,relative_DR_neg) :
		maxi.append(max(map(divide, a, b)))
	for a,b in zip(relative_ER,relative_ER_neg) :
		maxi.append(max(map(divide, a, b)))
	for a,b in zip(relative_IR,relative_IR_neg) :
		maxi.append(max(map(divide, a, b)))
		
	m = max(maxi)	
	ax1 = plt.subplot(gs[1])
	for a,b in zip(relative_DR,relative_DR_neg) :
		indexes1 = np.arange(Interdistance_maxValue + 1)
		#ax1 = fig.add_subplot(2,2,2)
		ax1.set_xlabel("base pairs between TGTCGG in DRs", size = 16)
		ax1.axis([-1, Interdistance_maxValue + 1, 0, m + 0.5])
		plt.plot(indexes1, map(divide, a, b), points[turn], alpha = 0.70, markersize = 13,markeredgewidth=0.0)
		plt.plot(indexes1, map(divide, a, b), lw=2)
		ax1.set_ylabel("DRn+ frequence / DRn- frequence", size = 16)
		plt.text(2, 3.5, "DRn frequence = DRn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$")
		l = plt.axhline(y = 1)
		for i in range (0,21) :
			plt.axvline(x = i, linewidth=0.1)
		turn = turn + 1
	turn = 0
	ax2 = plt.subplot(gs[2])
	for a, b,c in zip(relative_ER,relative_ER_neg,threshold)  :	
		indexes1 = range(Interdistance_maxValue + 1)
		#ax1 = fig.add_subplot(1,7,4)
		ax2.set_xlabel("base pairs between TGTCGG in ERs", size = 16)
		ax2.axis([-1, Interdistance_maxValue + 1, 0, m + 0.5])
		plt.plot(indexes1, map(divide, a, b), points[turn], label= r"threshold : "+str(c), alpha = 0.70, markersize = 13,markeredgewidth=0.0)
		plt.plot(indexes1, map(divide, a, b), lw=2)
		ax2.set_ylabel("ERn+ frequence / ERn- frequence", size = 16)
		plt.text(0, 3.5, "ERn frequence = ERn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$")
		l = plt.axhline(y = 1)
		plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc = "upper left", ncol = 2, mode="expand", borderaxespad=0.)
		turn = turn + 1
		for i in range (0,21) :
			plt.axvline(x = i, linewidth=0.1)
	turn = 0
	ax3 = plt.subplot(gs[3])
	for a,b in zip(relative_IR,relative_IR_neg) :
		indexes1 = range(Interdistance_maxValue + 1)
		#ax1 = fig.add_subplot(1,7,6)
		ax3.set_xlabel("base pairs between TGTCGG in IRs", size = 16)
		ax3.axis([-1, Interdistance_maxValue + 1, 0, m + 0.5])
		plt.plot(indexes1, map(divide, a, b), points[turn], alpha = 0.70, markersize = 13,markeredgewidth=0.0)
		plt.plot(indexes1, map(divide, a, b), lw=2)
		ax3.set_ylabel("IRn+ frequence / IRn- frequence", size = 16)
		plt.text(2, 3.5, "IRn frequence = IRn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$")
		l = plt.axhline(y = 1)
		turn = turn + 1
		for i in range (0,21) :
			plt.axvline(x = i, linewidth=0.1)
	plt.show()
	
	
#else :
	#for a in InterDR :
		#if factorTranscription == "LFY_matrix_19nucl" :
			#ax1 = fig.add_subplot(1,1,1)
		#else :	
			#ax1 = fig.add_subplot(1,3,1)
		#ax1.set_xlabel("base pairs between "+factorTranscription+" direct repeats", size = 16)
		#if not negative_sets and histo == True :
			#ax1.set_xlabel("base pairs between LFY 19 nucleotides matrix", size = 16)
			#indexes1 = range(20 + 1)
			#ax1.axis([0, 20 + 1, 0, 350])
			#ax1.bar(indexes1, a , width , color = 'cornflowerblue')
			#ax1.set_ylabel("occurences", color = 'cornflowerblue', size = 16)
			#indexes1 = np.arange(20 + 1)
			#plt.xticks(indexes1 + width * 0.5 , indexes1)
			#plt.text(8, 230, sys.argv)
		#else :
			#plt.plot(indexes1, a, lw=2)
			#ax1.axis([0, Interdistance_maxValue + 1, 0, 0.05])
			#ax1.set_ylabel("$DRn_+$ frequence", size = 16)	
	#if factorTranscription != "LFY_matrix_19nucl" :
		#for a, c in zip(InterER,threshold)  :
			#ax1 = fig.add_subplot(1,3,2)
			#ax1.set_xlabel("base pairs between "+factorTranscription+" everted repeats", size = 16)
			#if not negative_sets and histo == True :
				#indexes1 = range(Interdistance_maxValue + 1)
				#ax1.axis([0, Interdistance_maxValue + 1, 0, 350])
				#ax1.bar(indexes1, a , width , color = 'cornflowerblue')
				#indexes1 = np.arange(Interdistance_maxValue + 1)
				#plt.xticks(indexes1 + width * 0.5 , indexes1)
			#else :
				#plt.plot(indexes1, a, lw=2, label= "threshold : "+str(c))
				#ax1.axis([0, Interdistance_maxValue + 1, 0, 0.05])
				#ax1.set_ylabel("$ERn_+$ frequence", size = 16)
				#plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)	
	#if factorTranscription != "LFY_matrix_19nucl" :
		#for a in InterIR :	
			#ax1 = fig.add_subplot(1,3,3)
			#ax1.set_xlabel("base pairs between "+factorTranscription+" inverted repeats", size = 16)
			#if not negative_sets and histo == True :
				#indexes1 = range(Interdistance_maxValue + 1)
				#ax1.axis([0, Interdistance_maxValue + 1, 0, 350])
				#ax1.bar(indexes1, a , width , color = 'cornflowerblue')
				#indexes1 = np.arange(Interdistance_maxValue + 1)
				#plt.xticks(indexes1 + width * 0.5 , indexes1)
			#else:
				#plt.plot(indexes1, a, lw=2)
				#ax1.axis([0, Interdistance_maxValue + 1, 0, 0.05])
				#ax1.set_ylabel("$IRn_+$ frequence", size = 16)
	#plt.show()


