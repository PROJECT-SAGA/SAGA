
# -*- mode: python -*-

block_cipher = None

added_files = [ ('/home/ravel/SAGA/includes*', 'includes'), ]

a = Analysis(['SAGA.py'],
             pathex=['/home/ravel/SAGA/'],
             binaries=None,
             datas=added_files,
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='SAGA_linux',
          debug=False,
          strip=False,
          upx=True,
          console=False , icon='/home/ravel/SAGA/includes/icon.ico', resources=['/home/ravel/SAGA/includes/SAGA.ui'])


