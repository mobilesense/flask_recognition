import os

basedir = os.path.abspath(os.path.dirname(__file__))

#forms
WTF_CSRF_ENABLED = True
SECRET_KEY = 'you-will-never-guess'

SQLALCHEMY_DATABASE_URI = 'mysql://root:amine@localhost/yarilab_subscription_flask' 
SQLALCHEMY_COMMIT_ON_TEARDOWN = True

