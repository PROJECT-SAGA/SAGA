#!/usr/bin/python3.6
# -*- coding: utf-8 -*-
## @package SAGA.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################
import argparse
import sys, os

# Python QT modules
from PyQt5.QtWidgets import *
#from PyQt5.QtWidgets import QApplication, QFileSystemModel, QTableWidget, QFileDialog, QMessageBox, QTableWidgetItem
from PyQt5 import QtCore
from PyQt5.QtGui import QIcon, QPixmap
from PyQt5 import uic
from PyQt5.QtCore import Qt

#exit()
##################################################
## Variables Globales
version = '1.0'
VERSION_DATE = '09/04/2018'

##################################################
## Functions

def resource_path(relative):
	if hasattr(sys, "_MEIPASS"):
		return os.path.join(sys._MEIPASS, relative)
	return os.path.join(relative)


def existant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open by default.

	:param x: a file path
	:type x: str()
	:rtype: string
	:return: string

	"""
	if not os.path.exists(x):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))

	return x

filename = './includes/SAGA.ui'
myfontfile = resource_path(os.path.join(filename))
formSAGA,baseSAGA = uic.loadUiType(myfontfile)

def resource_path(relative_path):
	if hasattr(sys, '_MEIPASS'):
		return os.path.join(sys._MEIPASS, relative_path)
	return os.path.join(os.path.abspath("."), relative_path)

# ************************************* CLASSE SAGA Gestion de l'affichage et de la récupération de donnée ****************************************************
## @class SAGA
# @brief Classe principale, fenêtre principale
class SAGA( formSAGA, baseSAGA ):
	""" Classe principale qui est aussi la fenêtre principale de la GUI
	"""
	def __init__(self,app,parent=None):
		super(SAGA,self).__init__(parent)
		self.app = app
		self.ui = self
		self.ui.setupUi(self)
		self.ui.setFocus(True)
		self.ui.setWindowTitle("GUI_SAGA")
		self.ui.show()
		self.setFocusPolicy(QtCore.Qt.StrongFocus)
		self.activateWindow()

		self.initializationVariables()

		self.createWidgets()

	def initializationVariables(self):
		"""Initialize variables with defaults values"""
		self.dicoPath = {
						"vcf":None,
						"config":None,
						"reference":None,
						"trait":None,
						"SNPEFF":None,
						"GFF":None
					}
		self.dicoUINamePath = {
						"vcf":self.ui.vcfLineEdit,
						"config":self.ui.configLineEdit,
						"reference":self.ui.referenceLineEdit,
						"trait":self.ui.traitLineEdit,
						"SNPEFF":self.ui.SNPEFFLineEdit,
						"GFF":self.ui.GFFLineEdit
					}

		self.phasing = False

		self.SGEParam = ""

		# initialisation UI
		self.ui.runPushButton.setDisabled(True)
		self.ui.GFFFrame.hide()
		self.ui.SGEparamFrame.hide()

	def createWidgets(self):
		"""Mise en place du masquage des frames et des connections de l'interface"""
		self.setWindowIcon(QIcon(resource_path("./includes/icon.ico")))
		self.ui.statusbar.setStyleSheet("color: rgb(255, 107, 8);font: 8pt 'Arial Black';")

		## Edition des connect mandatory:
		self.ui.loadVCFPushButton.clicked.connect(lambda: self.loadFile(inputCat = "vcf", methode="clicked"))
		self.ui.vcfLineEdit.editingFinished.connect(lambda: self.loadFile(inputCat = "vcf", methode="write"))

		self.ui.loadConfigPushButton.clicked.connect(lambda: self.loadFile(inputCat = "config", methode="clicked"))
		self.ui.configLineEdit.editingFinished.connect(lambda: self.loadFile(inputCat = "config", methode="write"))

		self.ui.loadReferencePushButton.clicked.connect(lambda: self.loadFile(inputCat = "reference", methode="clicked"))
		self.ui.referenceLineEdit.editingFinished.connect(lambda: self.loadFile(inputCat = "reference", methode="write"))

		self.ui.loadTraitPushButton.clicked.connect(lambda: self.loadFile(inputCat = "trait", methode="clicked"))
		self.ui.traitLineEdit.editingFinished.connect(lambda: self.loadFile(inputCat = "trait", methode="write"))

		self.ui.runPushButton.clicked.connect(self.run)
		self.ui.resetPushButton.clicked.connect(self.reset)

		## Edition des options non mandatory:
		self.ui.gffsnpeffComboBox.currentIndexChanged[int].connect(self.actualizeGFFSNPEFFFrame)

		self.ui.loadSNPEFFPushButton.clicked.connect(lambda: self.loadFile(inputCat = "SNPEFF", methode="clicked"))
		self.ui.SNPEFFLineEdit.editingFinished.connect(lambda: self.loadFile(inputCat = "SNPEFF", methode="write"))

		self.ui.loadGFFPushButton.clicked.connect(lambda: self.loadFile(inputCat = "GFF", methode="clicked"))
		self.ui.GFFLineEdit.editingFinished.connect(lambda: self.loadFile(inputCat = "GFF", methode="write"))

		self.ui.phasingCheckBox.stateChanged[int].connect(self.actualizePhasing)

		self.ui.SGECheckBox.stateChanged[int].connect(self.actualizeSGEFrame)


	def actualizeGFFSNPEFFFrame(self):
		"""Change frame"""
		currentValue = str(self.ui.gffsnpeffComboBox.currentText())
		if "GFF" in currentValue:
			self.ui.GFFFrame.show()
			self.ui.SNPEFFFrame.hide()
			self.dicoPath["GFF"]= None
			self.dicoUINamePath["GFF"].setText("")
			self.actualizeRunButton()
		if "SNP" in currentValue:
			self.ui.GFFFrame.hide()
			self.ui.SNPEFFFrame.show()
			self.dicoPath["SNPEFF"]= None
			self.dicoUINamePath["SNPEFF"].setText("")
			self.actualizeRunButton()
		self.actualizeRunButton()


	def reset(self):
		"""To reset values"""
		for key in self.dicoUINamePath.keys():
			self.dicoPath[key]= None
			self.dicoUINamePath[inputCat].setText("")

		self.phasing = False
		self.ui.phasingCheckBox.setChecked(False)
		self.SGE = False
		self.ui.SGECheckBox.setChecked(False)
		self.SGEParam = ""

		self.actualizeRunButton()


	def actualizePhasing(self):
		"""change la valeur du choix quand changer"""
		if self.ui.phasingCheckBox.isChecked():
			self.phasing = True
		else:
			self.phasing = False
		self.actualizeRunButton()

	def actualizeSGEFrame(self):
		"""change la valeur du choix quand changer"""
		if self.ui.SGECheckBox.isChecked():
			self.ui.SGEparamFrame.show()
		else:
			self.ui.SGEparamFrame.hide()
			self.SGEParam = ""
			self.ui.SGEParamLineEdit.setText("")
		self.actualizeRunButton()


	def loadFile(self,inputCat = None, methode = None):
		"""Méthode qui permet de charger un fichier et afficher dans le plainText"""
		directoryToOpen = os.getcwd()
		if methode == "clicked":
			pathdir, _ = QFileDialog.getOpenFileName(self, caption="Load the file "+inputCat, directory=directoryToOpen)#filter="Text files (*.txt *.tab);;All (*.*)"
			self.dicoUINamePath[inputCat].setText(pathdir)
		elif methode == "write":
			pathdir = str(self.dicoUINamePath[inputCat].text())

		if pathdir != "" and os.path.isfile(pathdir):
			self.dicoPath[inputCat] = pathdir
		elif pathdir != "":
			self.dicoPath[inputCat]= None
			self.dicoUINamePath[inputCat].setText("")
			self.displayError(typeError="WARNING:",message = "\"%s\" is not a valid Path \n" % (pathdir))
		else:
			self.dicoPath[inputCat]= None
			self.dicoUINamePath[inputCat].setText("")
		self.actualizeRunButton()


	def actualizeRunButton(self):
		"""relache le bouton run si all mandatory"""
		if 	self.dicoPath["vcf"] != None and self.dicoPath["config"] != None and self.dicoPath["reference"] != None and self.dicoPath["trait"] != None and (self.dicoPath["SNPEFF"] != None or self.dicoPath["GFF"] != None ):
			self.ui.runPushButton.setEnabled(True)
		else:
			self.ui.runPushButton.setDisabled(True)


	def run(self):

		try:

			# for key, value in self.dicoPath.items():
				# print(key, value)
			# print(self.phasing)

			if self.phasing:
				phasingOpt = "-phasing"
			else:
				phasingOpt = ""

			if self.dicoPath["SNPEFF"] == None:
				database = "-gff {}".format(self.dicoPath["GFF"])
			elif self.dicoPath["GFF"] == None:
				database = "-snpeff {}".format(self.dicoPath["SNPEFF"])

			commandeLine = """
			saga.pl -vcf {}
			-conf {}
			-ref {}
			-trait {}
			{} {}""" .format(self.dicoPath["vcf"], self.dicoPath["config"], self.dicoPath["reference"], self.dicoPath["trait"], database, phasingOpt)



			if self.ui.SGECheckBox.isChecked():
				self.SGEParam = str(self.SGEParamLineEdit.text())
				commandeLine = "qsub {} {}".format(self.SGEParam, commandeLine)

			# print(commandeLine)

			self.displayError(typeError = "Command line:\n", message = commandeLine)


		except Exception as e:
				self.displayError(typeError = "ERROR:", message = "Error SAGA GUI....\n"+e)
				self.reset()


	def displayError(self,typeError,message):
		""" affiche les erreurs dans la zone de text"""
		if args.cmdMode:
			print(typeError,message)
		else:
			print(typeError,message)
			txtError = str(message)
			self.ui.statusbar.showMessage(txtError,7200)
			self.errorPopUp(typeError, message)

	def errorPopUp(self, typeError, message): ## Method to open a message box
		infoBox = QMessageBox()
		infoBox.setFixedSize(10000,100000)
		infoBox.setWindowTitle(typeError)
		if "ERROR" in typeError:
			infoBox.setIcon(QMessageBox.Critical)
		else:
			infoBox.setIcon(QMessageBox.Information)
		infoBox.setTextFormat(Qt.RichText)
		infoBox.setText("{} {}".format(typeError,message))
		infoBox.setTextInteractionFlags(Qt.TextSelectableByMouse)

		infoBox.exec_()


def main():

	#print sys.argv+["-cmd", "-gnome-terminal"]
	nargv = sys.argv

	# instanciation des objets principaux
	app = QApplication(nargv)
	myapp = SAGA(app)

	myapp.showMinimized()
	# les .app sous macos nécessitent cela pour que l'appli s'affiche en FG
	if "darwin" in sys.platform and ".app/" in sys.argv[0]:
		myapp.raise_()

	# lancement de la boucle Qt
	sys.exit(app.exec_())


# def cmd():
	# #print sys.argv+["-cmd", "-gnome-terminal"]
	# nargv = sys.argv
	# # instanciation des objets principaux
	# app = QApplication(nargv)
	# myapp = SAGA(app)

	# # Load info arguments to Class
	# # matrice file
	# myapp.matricePathFile = relativeToAbsolutePath(args.matriceParam)
	# # order file
	# myapp.orderPathFile = relativeToAbsolutePath(args.orderMatriceParam)
	# # PCA value
	# myapp.PCAvalue = args.pcaParam
	# # DA value
	# myapp.DAvalue = args.daParam
	# # pop min value
	# myapp.popMinValue = args.nbpopiParam
	# # pop max value
	# myapp.popMaxValue = args.nbpopmParam
	# # pgraph value
	# myapp.graphType = args.graphParam

	# # working dir path
	# workingDir = "/".join(relativeToAbsolutePath(args.orderMatriceParam).encode("utf-8").split("/")[:-1])+"/"
	# myapp.workingDir = workingDir.encode("utf-8")
	# # basename
	# basename = relativeToAbsolutePath(args.orderMatriceParam).encode("utf-8").split("/")[-1].split(".")[0]
	# myapp.basename = basename.encode("utf-8")
	# # pathFileOut dir path
	# pathFileOut = workingDir+basename+"/"
	# myapp.pathFileOut = pathFileOut.encode("utf-8")

	# # rm old folder
	# myapp.rmOld = args.rmOldParam


	# # Run programme
	# myapp.run()


if __name__ == '__main__':

	# Parameters recovery
	parser = argparse.ArgumentParser(prog='GUI_SAGA.py', description='''This Programme open GUI to produce SAGA script.\n
																				#If use on cluster you can run in commande line with option -c and args''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display GUI_SAGA.py version number and exit')
	parser.add_argument('-d', '--debug',action='store_true', help='enter verbose/debug mode', dest = "debug", default = "False")

	filesReq = parser.add_argument_group('Input mandatory infos for running if -c use')
	filesReq.add_argument('-c', '--cmd', action='store_true', dest = 'cmdMode', help = 'If used, programme run in CMD without interface')
	filesReq.add_argument('-m', '--mat', metavar="<filename>",type=existant_file, required=False, dest = 'matriceParam', help = 'matrice file path')
	filesReq.add_argument('-o', '--order', metavar="<filename>",type=existant_file, required=False, dest = 'orderMatriceParam', help = 'file with re-order name of matrice')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-pca', '--pcanum', metavar="<int>", required=False, default = "NULL", dest = 'pcaParam', help = 'Number value of PCA retains (default = NULL)')
	files.add_argument('-da', '--danum', metavar="<int>", required=False, default = "NULL", dest = 'daParam', help = 'Number value of DA retains (default = NULL)')
	files.add_argument('-pi', '--popi', metavar="<int>", type = int, default=2, required=False, dest = 'nbpopiParam', help = 'Number of pop Min (default = 2)')
	files.add_argument('-pm', '--popm', metavar="<int>", type = int, default=10, required=False, dest = 'nbpopmParam', help = 'Number of pop Max (default = 10)')	## Check parameters
	files.add_argument('-r', '--rm', metavar="<True/False>", choices=("True","False"), default="False", required=False, dest = 'rmOldParam', help = 'if directory exist remove (default = False)')	## Check parameters
	files.add_argument('-g', '--graph', metavar="<1/2/3>", choices=("1","2","3"), default="1", required=False, dest = 'graphParam', help = 'type of graph (default = 1)')	## Check parameters
	args = parser.parse_args()

	if args.cmdMode:
		if args.matriceParam == "" or args.orderMatriceParam == "":
			print("ERROR: You must enter require arguments")
			exit()
		cmd()
	else:
		# run interface
		main()
