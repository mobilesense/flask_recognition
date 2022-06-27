import os

basedir = os.path.abspath(os.path.dirname(__file__))

#forms
WTF_CSRF_ENABLED = True
SECRET_KEY = 'you-will-never-guess'

SQLALCHEMY_DATABASE_URI = 'mysql://root:amine@localhost/yarilab_flask' 
SQLALCHEMY_COMMIT_ON_TEARDOWN = True

#mail
MAIL_SERVER = 'smtp.yarilab.com'
MAIL_PORT = 25
MAIL_USE_SSL = True
MAIL_USERNAME = 'amine'
MAIL_PASSWORD = 'password'

#pagination
PER_PAGE = 5
LINK_SIZE = 'sm'
# decide whether or not a single page returns pagination
SHOW_SINGLE_PAGE = False
CSS_FRAMEWORK = 'bootstrap3'

#search
MAX_SEARCH_RESULTS = 50

#data per user 
DATA_FOLDER = './data'
APP_BASE_DIR = '/home/amine/Documents/src/flask_recognition/'

#google login 
GOOGLE_LOGIN_CLIENT_ID='933803569112-3qri3m1p5a9ghv6c37q03vqh9ibnd54c.apps.googleusercontent.com ',
GOOGLE_LOGIN_CLIENT_SECRET='EHDScAcuP00UTIn4zQkHzfSC',
GOOGLE_LOGIN_REDIRECT_URI='https://www.yarilab.com/oauth2callback'#'https://www.yarilab.com/oauth2callback'

# Excess toleration in months
MAX_TOLERATE_EXCESS = 1 #within this period the API is running.


