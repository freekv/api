# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'org/helioviewer/InstallWizard.ui'
#
# Created: Tue May  3 05:58:08 2011
#      by: PyQt4 UI code generator 4.8.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_InstallWizard(object):
    def setupUi(self, InstallWizard):
        InstallWizard.setObjectName(_fromUtf8("InstallWizard"))
        InstallWizard.resize(539, 391)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(_fromUtf8(":/Logos/hv.ico")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        InstallWizard.setWindowIcon(icon)
        InstallWizard.setWizardStyle(QtGui.QWizard.AeroStyle)
        InstallWizard.setOptions(QtGui.QWizard.DisabledBackButtonOnLastPage|QtGui.QWizard.IndependentPages)
        self.introPage = QtGui.QWizardPage()
        self.introPage.setSubTitle(_fromUtf8(""))
        self.introPage.setObjectName(_fromUtf8("introPage"))
        self.verticalLayoutWidget = QtGui.QWidget(self.introPage)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(10, 0, 511, 141))
        self.verticalLayoutWidget.setObjectName(_fromUtf8("verticalLayoutWidget"))
        self.introVBox = QtGui.QVBoxLayout(self.verticalLayoutWidget)
        self.introVBox.setMargin(0)
        self.introVBox.setObjectName(_fromUtf8("introVBox"))
        self.greetingMsg = QtGui.QLabel(self.verticalLayoutWidget)
        self.greetingMsg.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.greetingMsg.setWordWrap(True)
        self.greetingMsg.setIndent(13)
        self.greetingMsg.setOpenExternalLinks(True)
        self.greetingMsg.setObjectName(_fromUtf8("greetingMsg"))
        self.introVBox.addWidget(self.greetingMsg)
        InstallWizard.addPage(self.introPage)
        self.dbAdminPage = QtGui.QWizardPage()
        self.dbAdminPage.setObjectName(_fromUtf8("dbAdminPage"))
        self.verticalLayoutWidget_4 = QtGui.QWidget(self.dbAdminPage)
        self.verticalLayoutWidget_4.setGeometry(QtCore.QRect(0, 0, 521, 288))
        self.verticalLayoutWidget_4.setObjectName(_fromUtf8("verticalLayoutWidget_4"))
        self.dbAdminPageOuter = QtGui.QVBoxLayout(self.verticalLayoutWidget_4)
        self.dbAdminPageOuter.setMargin(0)
        self.dbAdminPageOuter.setObjectName(_fromUtf8("dbAdminPageOuter"))
        self.dbAdminPageDesc = QtGui.QLabel(self.verticalLayoutWidget_4)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.dbAdminPageDesc.sizePolicy().hasHeightForWidth())
        self.dbAdminPageDesc.setSizePolicy(sizePolicy)
        self.dbAdminPageDesc.setMinimumSize(QtCore.QSize(0, 70))
        self.dbAdminPageDesc.setWordWrap(True)
        self.dbAdminPageDesc.setObjectName(_fromUtf8("dbAdminPageDesc"))
        self.dbAdminPageOuter.addWidget(self.dbAdminPageDesc)
        self.dbAdminForm = QtGui.QFormLayout()
        self.dbAdminForm.setFieldGrowthPolicy(QtGui.QFormLayout.AllNonFixedFieldsGrow)
        self.dbAdminForm.setObjectName(_fromUtf8("dbAdminForm"))
        self.dbTypeLbl = QtGui.QLabel(self.verticalLayoutWidget_4)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.dbTypeLbl.sizePolicy().hasHeightForWidth())
        self.dbTypeLbl.setSizePolicy(sizePolicy)
        self.dbTypeLbl.setMinimumSize(QtCore.QSize(0, 57))
        self.dbTypeLbl.setAlignment(QtCore.Qt.AlignBottom|QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft)
        self.dbTypeLbl.setIndent(20)
        self.dbTypeLbl.setObjectName(_fromUtf8("dbTypeLbl"))
        self.dbAdminForm.setWidget(0, QtGui.QFormLayout.LabelRole, self.dbTypeLbl)
        self.dbType = QtGui.QWidget(self.verticalLayoutWidget_4)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.dbType.sizePolicy().hasHeightForWidth())
        self.dbType.setSizePolicy(sizePolicy)
        self.dbType.setMinimumSize(QtCore.QSize(0, 50))
        self.dbType.setBaseSize(QtCore.QSize(0, 0))
        self.dbType.setObjectName(_fromUtf8("dbType"))
        self.postgresRadioBtn = QtGui.QRadioButton(self.dbType)
        self.postgresRadioBtn.setGeometry(QtCore.QRect(140, 10, 392, 41))
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.postgresRadioBtn.sizePolicy().hasHeightForWidth())
        self.postgresRadioBtn.setSizePolicy(sizePolicy)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(_fromUtf8(":/Databases/postgres_logo-elephant-64.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.postgresRadioBtn.setIcon(icon1)
        self.postgresRadioBtn.setIconSize(QtCore.QSize(32, 32))
        self.postgresRadioBtn.setObjectName(_fromUtf8("postgresRadioBtn"))
        self.mysqlRadioBtn = QtGui.QRadioButton(self.dbType)
        self.mysqlRadioBtn.setGeometry(QtCore.QRect(0, 10, 91, 37))
        self.mysqlRadioBtn.setText(_fromUtf8(""))
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(_fromUtf8(":/Databases/mysql_logo.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.mysqlRadioBtn.setIcon(icon2)
        self.mysqlRadioBtn.setIconSize(QtCore.QSize(63, 32))
        self.mysqlRadioBtn.setChecked(True)
        self.mysqlRadioBtn.setObjectName(_fromUtf8("mysqlRadioBtn"))
        self.dbAdminForm.setWidget(0, QtGui.QFormLayout.FieldRole, self.dbType)
        self.dbAdminUserNameLbl = QtGui.QLabel(self.verticalLayoutWidget_4)
        self.dbAdminUserNameLbl.setIndent(20)
        self.dbAdminUserNameLbl.setObjectName(_fromUtf8("dbAdminUserNameLbl"))
        self.dbAdminForm.setWidget(1, QtGui.QFormLayout.LabelRole, self.dbAdminUserNameLbl)
        self.dbAdminUserName = QtGui.QLineEdit(self.verticalLayoutWidget_4)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.dbAdminUserName.sizePolicy().hasHeightForWidth())
        self.dbAdminUserName.setSizePolicy(sizePolicy)
        self.dbAdminUserName.setObjectName(_fromUtf8("dbAdminUserName"))
        self.dbAdminForm.setWidget(1, QtGui.QFormLayout.FieldRole, self.dbAdminUserName)
        spacerItem = QtGui.QSpacerItem(20, 10, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.dbAdminForm.setItem(2, QtGui.QFormLayout.FieldRole, spacerItem)
        self.dbAdminPasswordLbl = QtGui.QLabel(self.verticalLayoutWidget_4)
        self.dbAdminPasswordLbl.setIndent(20)
        self.dbAdminPasswordLbl.setObjectName(_fromUtf8("dbAdminPasswordLbl"))
        self.dbAdminForm.setWidget(3, QtGui.QFormLayout.LabelRole, self.dbAdminPasswordLbl)
        self.dbAdminPassword = QtGui.QLineEdit(self.verticalLayoutWidget_4)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.dbAdminPassword.sizePolicy().hasHeightForWidth())
        self.dbAdminPassword.setSizePolicy(sizePolicy)
        self.dbAdminPassword.setEchoMode(QtGui.QLineEdit.Password)
        self.dbAdminPassword.setObjectName(_fromUtf8("dbAdminPassword"))
        self.dbAdminForm.setWidget(3, QtGui.QFormLayout.FieldRole, self.dbAdminPassword)
        self.dbAdminStatus = QtGui.QLabel(self.verticalLayoutWidget_4)
        self.dbAdminStatus.setMinimumSize(QtCore.QSize(0, 50))
        self.dbAdminStatus.setText(_fromUtf8(""))
        self.dbAdminStatus.setWordWrap(True)
        self.dbAdminStatus.setObjectName(_fromUtf8("dbAdminStatus"))
        self.dbAdminForm.setWidget(4, QtGui.QFormLayout.FieldRole, self.dbAdminStatus)
        self.dbAdminPageOuter.addLayout(self.dbAdminForm)
        InstallWizard.addPage(self.dbAdminPage)
        self.hvDatabaseSetupPage = QtGui.QWizardPage()
        self.hvDatabaseSetupPage.setObjectName(_fromUtf8("hvDatabaseSetupPage"))
        self.verticalLayoutWidget_3 = QtGui.QWidget(self.hvDatabaseSetupPage)
        self.verticalLayoutWidget_3.setGeometry(QtCore.QRect(0, 0, 521, 271))
        self.verticalLayoutWidget_3.setObjectName(_fromUtf8("verticalLayoutWidget_3"))
        self.hvDatabaseSetupPageOuter = QtGui.QVBoxLayout(self.verticalLayoutWidget_3)
        self.hvDatabaseSetupPageOuter.setMargin(0)
        self.hvDatabaseSetupPageOuter.setObjectName(_fromUtf8("hvDatabaseSetupPageOuter"))
        self.hvDatabaseSetupPageDesc = QtGui.QLabel(self.verticalLayoutWidget_3)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.hvDatabaseSetupPageDesc.sizePolicy().hasHeightForWidth())
        self.hvDatabaseSetupPageDesc.setSizePolicy(sizePolicy)
        self.hvDatabaseSetupPageDesc.setMinimumSize(QtCore.QSize(0, 70))
        self.hvDatabaseSetupPageDesc.setWordWrap(True)
        self.hvDatabaseSetupPageDesc.setObjectName(_fromUtf8("hvDatabaseSetupPageDesc"))
        self.hvDatabaseSetupPageOuter.addWidget(self.hvDatabaseSetupPageDesc)
        self.hvDatabaseSetupForm = QtGui.QFormLayout()
        self.hvDatabaseSetupForm.setFieldGrowthPolicy(QtGui.QFormLayout.AllNonFixedFieldsGrow)
        self.hvDatabaseSetupForm.setObjectName(_fromUtf8("hvDatabaseSetupForm"))
        self.hvUserNameLbl = QtGui.QLabel(self.verticalLayoutWidget_3)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.hvUserNameLbl.sizePolicy().hasHeightForWidth())
        self.hvUserNameLbl.setSizePolicy(sizePolicy)
        self.hvUserNameLbl.setIndent(20)
        self.hvUserNameLbl.setObjectName(_fromUtf8("hvUserNameLbl"))
        self.hvDatabaseSetupForm.setWidget(2, QtGui.QFormLayout.LabelRole, self.hvUserNameLbl)
        self.hvUserName = QtGui.QLineEdit(self.verticalLayoutWidget_3)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.hvUserName.sizePolicy().hasHeightForWidth())
        self.hvUserName.setSizePolicy(sizePolicy)
        self.hvUserName.setText(_fromUtf8(""))
        self.hvUserName.setObjectName(_fromUtf8("hvUserName"))
        self.hvDatabaseSetupForm.setWidget(2, QtGui.QFormLayout.FieldRole, self.hvUserName)
        self.hvPasswordLbl = QtGui.QLabel(self.verticalLayoutWidget_3)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.hvPasswordLbl.sizePolicy().hasHeightForWidth())
        self.hvPasswordLbl.setSizePolicy(sizePolicy)
        self.hvPasswordLbl.setIndent(20)
        self.hvPasswordLbl.setObjectName(_fromUtf8("hvPasswordLbl"))
        self.hvDatabaseSetupForm.setWidget(4, QtGui.QFormLayout.LabelRole, self.hvPasswordLbl)
        self.hvPassword = QtGui.QLineEdit(self.verticalLayoutWidget_3)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.hvPassword.sizePolicy().hasHeightForWidth())
        self.hvPassword.setSizePolicy(sizePolicy)
        self.hvPassword.setText(_fromUtf8(""))
        self.hvPassword.setEchoMode(QtGui.QLineEdit.Normal)
        self.hvPassword.setObjectName(_fromUtf8("hvPassword"))
        self.hvDatabaseSetupForm.setWidget(4, QtGui.QFormLayout.FieldRole, self.hvPassword)
        spacerItem1 = QtGui.QSpacerItem(20, 10, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.hvDatabaseSetupForm.setItem(3, QtGui.QFormLayout.FieldRole, spacerItem1)
        self.hvDatabaseSetupStatus = QtGui.QLabel(self.verticalLayoutWidget_3)
        self.hvDatabaseSetupStatus.setText(_fromUtf8(""))
        self.hvDatabaseSetupStatus.setWordWrap(True)
        self.hvDatabaseSetupStatus.setObjectName(_fromUtf8("hvDatabaseSetupStatus"))
        self.hvDatabaseSetupForm.setWidget(5, QtGui.QFormLayout.FieldRole, self.hvDatabaseSetupStatus)
        self.hvDatabaseNameLbl = QtGui.QLabel(self.verticalLayoutWidget_3)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.hvDatabaseNameLbl.sizePolicy().hasHeightForWidth())
        self.hvDatabaseNameLbl.setSizePolicy(sizePolicy)
        self.hvDatabaseNameLbl.setIndent(20)
        self.hvDatabaseNameLbl.setObjectName(_fromUtf8("hvDatabaseNameLbl"))
        self.hvDatabaseSetupForm.setWidget(0, QtGui.QFormLayout.LabelRole, self.hvDatabaseNameLbl)
        self.hvDatabaseName = QtGui.QLineEdit(self.verticalLayoutWidget_3)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.hvDatabaseName.sizePolicy().hasHeightForWidth())
        self.hvDatabaseName.setSizePolicy(sizePolicy)
        self.hvDatabaseName.setText(_fromUtf8(""))
        self.hvDatabaseName.setObjectName(_fromUtf8("hvDatabaseName"))
        self.hvDatabaseSetupForm.setWidget(0, QtGui.QFormLayout.FieldRole, self.hvDatabaseName)
        spacerItem2 = QtGui.QSpacerItem(20, 10, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.hvDatabaseSetupForm.setItem(1, QtGui.QFormLayout.FieldRole, spacerItem2)
        self.hvDatabaseSetupPageOuter.addLayout(self.hvDatabaseSetupForm)
        InstallWizard.addPage(self.hvDatabaseSetupPage)
        self.jp2ArchivePage = QtGui.QWizardPage()
        self.jp2ArchivePage.setObjectName(_fromUtf8("jp2ArchivePage"))
        self.verticalLayoutWidget_5 = QtGui.QWidget(self.jp2ArchivePage)
        self.verticalLayoutWidget_5.setGeometry(QtCore.QRect(-1, -1, 521, 271))
        self.verticalLayoutWidget_5.setObjectName(_fromUtf8("verticalLayoutWidget_5"))
        self.jp2ArchivePageOuter = QtGui.QVBoxLayout(self.verticalLayoutWidget_5)
        self.jp2ArchivePageOuter.setMargin(0)
        self.jp2ArchivePageOuter.setObjectName(_fromUtf8("jp2ArchivePageOuter"))
        self.jp2ArchivePageDesc = QtGui.QLabel(self.verticalLayoutWidget_5)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.jp2ArchivePageDesc.sizePolicy().hasHeightForWidth())
        self.jp2ArchivePageDesc.setSizePolicy(sizePolicy)
        self.jp2ArchivePageDesc.setMinimumSize(QtCore.QSize(0, 70))
        self.jp2ArchivePageDesc.setObjectName(_fromUtf8("jp2ArchivePageDesc"))
        self.jp2ArchivePageOuter.addWidget(self.jp2ArchivePageDesc)
        self.jp2ArchivePageHBox = QtGui.QHBoxLayout()
        self.jp2ArchivePageHBox.setContentsMargins(-1, -1, 15, -1)
        self.jp2ArchivePageHBox.setObjectName(_fromUtf8("jp2ArchivePageHBox"))
        self.jp2RootDirLbl = QtGui.QLabel(self.verticalLayoutWidget_5)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.jp2RootDirLbl.sizePolicy().hasHeightForWidth())
        self.jp2RootDirLbl.setSizePolicy(sizePolicy)
        self.jp2RootDirLbl.setMaximumSize(QtCore.QSize(80, 16777215))
        self.jp2RootDirLbl.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.jp2RootDirLbl.setWordWrap(True)
        self.jp2RootDirLbl.setObjectName(_fromUtf8("jp2RootDirLbl"))
        self.jp2ArchivePageHBox.addWidget(self.jp2RootDirLbl)
        self.jp2RootDirInput = QtGui.QLineEdit(self.verticalLayoutWidget_5)
        self.jp2RootDirInput.setObjectName(_fromUtf8("jp2RootDirInput"))
        self.jp2ArchivePageHBox.addWidget(self.jp2RootDirInput)
        self.jp2BrowseBtn = QtGui.QPushButton(self.verticalLayoutWidget_5)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.jp2BrowseBtn.sizePolicy().hasHeightForWidth())
        self.jp2BrowseBtn.setSizePolicy(sizePolicy)
        self.jp2BrowseBtn.setObjectName(_fromUtf8("jp2BrowseBtn"))
        self.jp2ArchivePageHBox.addWidget(self.jp2BrowseBtn)
        self.jp2ArchivePageOuter.addLayout(self.jp2ArchivePageHBox)
        self.jp2ArchiveStatus = QtGui.QLabel(self.verticalLayoutWidget_5)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.jp2ArchiveStatus.sizePolicy().hasHeightForWidth())
        self.jp2ArchiveStatus.setSizePolicy(sizePolicy)
        self.jp2ArchiveStatus.setMinimumSize(QtCore.QSize(0, 50))
        self.jp2ArchiveStatus.setText(_fromUtf8(""))
        self.jp2ArchiveStatus.setWordWrap(True)
        self.jp2ArchiveStatus.setObjectName(_fromUtf8("jp2ArchiveStatus"))
        self.jp2ArchivePageOuter.addWidget(self.jp2ArchiveStatus)
        InstallWizard.addPage(self.jp2ArchivePage)
        self.installingPage = QtGui.QWizardPage()
        self.installingPage.setObjectName(_fromUtf8("installingPage"))
        self.verticalLayoutWidget_2 = QtGui.QWidget(self.installingPage)
        self.verticalLayoutWidget_2.setGeometry(QtCore.QRect(0, 0, 521, 241))
        self.verticalLayoutWidget_2.setObjectName(_fromUtf8("verticalLayoutWidget_2"))
        self.installingPageOuter = QtGui.QVBoxLayout(self.verticalLayoutWidget_2)
        self.installingPageOuter.setContentsMargins(10, -1, 20, -1)
        self.installingPageOuter.setObjectName(_fromUtf8("installingPageOuter"))
        self.label = QtGui.QLabel(self.verticalLayoutWidget_2)
        self.label.setAlignment(QtCore.Qt.AlignBottom|QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft)
        self.label.setIndent(0)
        self.label.setObjectName(_fromUtf8("label"))
        self.installingPageOuter.addWidget(self.label)
        self.statusMsg = QtGui.QLabel(self.verticalLayoutWidget_2)
        self.statusMsg.setText(_fromUtf8(""))
        self.statusMsg.setIndent(15)
        self.statusMsg.setObjectName(_fromUtf8("statusMsg"))
        self.installingPageOuter.addWidget(self.statusMsg)
        self.installProgress = QtGui.QProgressBar(self.verticalLayoutWidget_2)
        self.installProgress.setEnabled(True)
        self.installProgress.setProperty(_fromUtf8("value"), 0)
        self.installProgress.setObjectName(_fromUtf8("installProgress"))
        self.installingPageOuter.addWidget(self.installProgress)
        self.widget = QtGui.QWidget(self.verticalLayoutWidget_2)
        self.widget.setMaximumSize(QtCore.QSize(16777215, 50))
        self.widget.setObjectName(_fromUtf8("widget"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.widget)
        self.horizontalLayout.setMargin(0)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        spacerItem3 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem3)
        self.startProcessingBtn = QtGui.QPushButton(self.widget)
        self.startProcessingBtn.setObjectName(_fromUtf8("startProcessingBtn"))
        self.horizontalLayout.addWidget(self.startProcessingBtn)
        self.installingPageOuter.addWidget(self.widget)
        InstallWizard.addPage(self.installingPage)
        self.finishedPage = QtGui.QWizardPage()
        self.finishedPage.setObjectName(_fromUtf8("finishedPage"))
        self.finishedPageDesc = QtGui.QLabel(self.finishedPage)
        self.finishedPageDesc.setGeometry(QtCore.QRect(20, 0, 501, 51))
        self.finishedPageDesc.setObjectName(_fromUtf8("finishedPageDesc"))
        InstallWizard.addPage(self.finishedPage)

        self.retranslateUi(InstallWizard)
        QtCore.QMetaObject.connectSlotsByName(InstallWizard)
        InstallWizard.setTabOrder(self.mysqlRadioBtn, self.postgresRadioBtn)
        InstallWizard.setTabOrder(self.postgresRadioBtn, self.dbAdminUserName)
        InstallWizard.setTabOrder(self.dbAdminUserName, self.dbAdminPassword)
        InstallWizard.setTabOrder(self.dbAdminPassword, self.hvDatabaseName)
        InstallWizard.setTabOrder(self.hvDatabaseName, self.hvUserName)
        InstallWizard.setTabOrder(self.hvUserName, self.hvPassword)
        InstallWizard.setTabOrder(self.hvPassword, self.jp2RootDirInput)
        InstallWizard.setTabOrder(self.jp2RootDirInput, self.jp2BrowseBtn)
        InstallWizard.setTabOrder(self.jp2BrowseBtn, self.startProcessingBtn)

    def retranslateUi(self, InstallWizard):
        InstallWizard.setWindowTitle(QtGui.QApplication.translate("InstallWizard", "Helioviewer Database Installation", None, QtGui.QApplication.UnicodeUTF8))
        self.introPage.setTitle(QtGui.QApplication.translate("InstallWizard", "Introduction", None, QtGui.QApplication.UnicodeUTF8))
        self.greetingMsg.setText(QtGui.QApplication.translate("InstallWizard", "This wizard will help setup the required database structure used by the Helioviewer project. We assume that you have an existing database installation (either MySQL or PostgreSQL), and that you have the proper access priveges in order to create new schemas and users in that system. Click \"Next\" to begin!", None, QtGui.QApplication.UnicodeUTF8))
        self.dbAdminPage.setTitle(QtGui.QApplication.translate("InstallWizard", "Database Setup", None, QtGui.QApplication.UnicodeUTF8))
        self.dbAdminPage.setSubTitle(QtGui.QApplication.translate("InstallWizard", "Existing database information", None, QtGui.QApplication.UnicodeUTF8))
        self.dbAdminPageDesc.setText(QtGui.QApplication.translate("InstallWizard", "Select the type of database you wish to use and the login information for an administrative user in that system", None, QtGui.QApplication.UnicodeUTF8))
        self.dbTypeLbl.setText(QtGui.QApplication.translate("InstallWizard", "Database type:", None, QtGui.QApplication.UnicodeUTF8))
        self.postgresRadioBtn.setText(QtGui.QApplication.translate("InstallWizard", "PostgreSQL", None, QtGui.QApplication.UnicodeUTF8))
        self.dbAdminUserNameLbl.setText(QtGui.QApplication.translate("InstallWizard", "Username:", None, QtGui.QApplication.UnicodeUTF8))
        self.dbAdminPasswordLbl.setText(QtGui.QApplication.translate("InstallWizard", "Password:", None, QtGui.QApplication.UnicodeUTF8))
        self.hvDatabaseSetupPage.setTitle(QtGui.QApplication.translate("InstallWizard", "Database Setup", None, QtGui.QApplication.UnicodeUTF8))
        self.hvDatabaseSetupPage.setSubTitle(QtGui.QApplication.translate("InstallWizard", "New database information", None, QtGui.QApplication.UnicodeUTF8))
        self.hvDatabaseSetupPageDesc.setText(QtGui.QApplication.translate("InstallWizard", "Select a name for the new database to be created, and the username and password you wish to use to access this database.", None, QtGui.QApplication.UnicodeUTF8))
        self.hvUserNameLbl.setText(QtGui.QApplication.translate("InstallWizard", "Username:", None, QtGui.QApplication.UnicodeUTF8))
        self.hvPasswordLbl.setText(QtGui.QApplication.translate("InstallWizard", "Password:", None, QtGui.QApplication.UnicodeUTF8))
        self.hvDatabaseNameLbl.setText(QtGui.QApplication.translate("InstallWizard", "Database:", None, QtGui.QApplication.UnicodeUTF8))
        self.jp2ArchivePage.setTitle(QtGui.QApplication.translate("InstallWizard", "Image Archive", None, QtGui.QApplication.UnicodeUTF8))
        self.jp2ArchivePage.setSubTitle(QtGui.QApplication.translate("InstallWizard", "JPEG 2000 Image Archive Location", None, QtGui.QApplication.UnicodeUTF8))
        self.jp2ArchivePageDesc.setText(QtGui.QApplication.translate("InstallWizard", "Enter the location of the JPEG 2000 image archive root directory:", None, QtGui.QApplication.UnicodeUTF8))
        self.jp2RootDirLbl.setText(QtGui.QApplication.translate("InstallWizard", "JP2 Archive Location:", None, QtGui.QApplication.UnicodeUTF8))
        self.jp2RootDirInput.setText(QtGui.QApplication.translate("InstallWizard", "/var/www/jp2", None, QtGui.QApplication.UnicodeUTF8))
        self.jp2BrowseBtn.setText(QtGui.QApplication.translate("InstallWizard", "Browse", None, QtGui.QApplication.UnicodeUTF8))
        self.installingPage.setTitle(QtGui.QApplication.translate("InstallWizard", "Installing Helioviewer Database Schema", None, QtGui.QApplication.UnicodeUTF8))
        self.installingPage.setSubTitle(QtGui.QApplication.translate("InstallWizard", "Processing JPEG 2000 Archive", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("InstallWizard", "<b>Status:</b>", None, QtGui.QApplication.UnicodeUTF8))
        self.startProcessingBtn.setText(QtGui.QApplication.translate("InstallWizard", "Start", None, QtGui.QApplication.UnicodeUTF8))
        self.finishedPage.setTitle(QtGui.QApplication.translate("InstallWizard", "Finished!", None, QtGui.QApplication.UnicodeUTF8))
        self.finishedPageDesc.setText(QtGui.QApplication.translate("InstallWizard", "You have successfully installed the Helioviewer database.", None, QtGui.QApplication.UnicodeUTF8))

from .resources.hv_rc import *