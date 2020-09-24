import re
import os
from os import path
import subprocess


potentialString= r'potential(psi_real,psi_imag,x){ psi_real*psi_real + psi_imag*psi_imag + 0.5 *( x[0]**2 + x[1]**2)}'

def parseLocalPotential( code , coordinate="cylindrical"):
	'''
	Returns fragment of C++ code inside the evaluation loop
	for function evaluation
	'''
	if coordinate is not "cylindrical":
		raise NotImplementedError("Coordinates " + coordinate )
	pattern=r'\s*(\w+)\s*\(\s*(\w+)\s*,\s*(\w+)\s*,\s*(\w+)\s*\)\s*\{(.*)\}'
	m=re.match(pattern,code)
	name,phi_real,phi_imag,x,body=m.groups()
	#psi_real = m.goup(2)
	#psi_imag = m.group(3)
	#x = m.group(4)
	#body = m.group(5)
	
	return expandVariables(body,phi_real,phi_imag,x)

def expandVariables(body,psi_real,psi_imag,x):

	body=re.sub(psi_real, r"phi_old_real(i,j,0)",body)
	body=re.sub(psi_imag, r"phi_old_imag(i,j,0)",body)
	body=re.sub(x + r"\[0\]","r",body)
	body=re.sub(x + r"\[1\]","z",body)

	return createAssignments(body)

def createAssignments(body):
	code =r"Real tmp = " + body + ";\n"
	code+= r"phi_new_real(i,j,0)+= tmp*phi_old_real(i,j,0);" + "\n"
	code+= r"phi_new_imag(i,j,0)+= tmp*phi_old_imag(i,j,0);" + "\n"
	return code

def compileLocalPotential(code,coordinate="cylindrical"):
	code=parseLocalPotential(code,coordinate)

	code_dir=path.dirname( __file__ ) + "/.."
	template_filename=code_dir + "/evaluate_templates/evaluate_cylindrical_template.cpp"
	script_file=code_dir + "/set_evaluate.sh"
	
	with open(template_filename) as f:	
		out=f.read()
	out=re.sub("LOCAL_ASSIGNMENT_MARKER",code,out)

	with open(code_dir + "/evaluate.cpp","w") as f:
		f.write(out)

	compile_output=subprocess.check_output(["bash" ,"-c","make install -C " + code_dir + "/build"])

	return compile_output




