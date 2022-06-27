import flask
from flask import url_for, request, abort, jsonify, g, session, flash, redirect
from flask.ext.login import login_user, logout_user, current_user, \
    login_required
from itsdangerous import URLSafeSerializer, BadSignature    
from flask.ext.mail import Message    
from flask.ext.paginate import Pagination  
from sqlalchemy.sql import func
from sqlalchemy import desc
from sqlalchemy import and_
from app import app, db, login_manager, mail #, googlelogin
from .models import User, Project, Visual, Query, Mobile, Usage, Billing, Notification
from .forms import LoginForm, RegistrationForm, ResetForm, \
    EditUserInfoForm, ChangePassForm, NewProjectForm, EditProjectForm,\
    EditVisualForm, SearchForm, EditMobileForm
from app.globals import mail_msg_body, plans_limit
from pytz import timezone
from werkzeug import secure_filename
#yarilab core api
#import core

#divers
from calendar import monthrange
import datetime
from dateutil.relativedelta import relativedelta
import hashlib
import logging
import os
import shutil
import sys
import re
import pdb

"""  
remark: 
- login_form and reset_form are shared for all "non login required pages" because they are in header.
- do not pass login_form and reset_form to "login_required pages" -> no effect
"""

ALLOWED_EXTENSIONS = set(['png', 'jpg', 'jpeg'])

def get_notifications(limit=5):
    if g.user.is_authenticated():
        if limit != "all":
            return Notification.query.filter_by(user_id = current_user.id).order_by(desc(Notification.creation_date)).limit(5).all()
        else:
            return Notification.query.filter_by(user_id = current_user.id).order_by(desc(Notification.creation_date)).offset(5).all()        
    else:
        return []    

@login_manager.user_loader
#@googlelogin.user_loader
def load_user(userid):
    return User.query.get(userid)

@app.before_request
def before_request():
    g.user = current_user

def root_dir():  
    return os.path.abspath(os.path.dirname(__file__))

def google_login_urls():
    #urls = [googlelogin.login_url(approval_prompt='force'), googlelogin.login_url(approval_prompt='force',params=dict(extra='large-fries')), googlelogin.login_url(approval_prompt='force',scopes=['https://www.googleapis.com/auth/drive'],access_type='offline')]
    urls = []
    return urls
        
@app.route('/search/project', methods=['POST'])
@login_required
def project_search():
    project_search_form = SearchForm(request.form)
    if not project_search_form.validate_on_submit():
        return redirect(url_for('index'))
    return redirect(url_for('manage_database', query=project_search_form.search.data))

@app.route('/oauth2callback')
#@googlelogin.oauth2callback
def google_login(token, userinfo, **params):
    user = users[userinfo['id']] = User(userinfo)
    login_user(user)
    session['token'] = json.dumps(token)
    session['extra'] = params.get('extra')
    return redirect(request.values.get('next'))
            
@app.route("/login", methods=["POST"])
def login():
    if g.user is not None and g.user.is_authenticated():    
        return redirect(url_for("index"))
    form = LoginForm(request.form)        
    if form.validate_on_submit():
        user = User.query.filter_by(username=form.username.data).first()
        if user and user.verify_password(form.password.data):                    
            if user.active == True:                
                session['remember_me'] = form.remember_me.data     
                session.permanent = True
                login_user(user) 
                """user.last_session = datetime.datetime.now(timezone('UTC'))  
                db.session.commit()  """
                return redirect(request.values.get('next'))
            else:
                flash("Check your mail box for account activation link.", 'login_error')                
                return redirect(request.values.get('next'))
        else:           
            flash("Please enter a correct username and password.", 'login_error')            
            return redirect(request.values.get('next')) #, login_form=form
    if len(form.username.errors)>0: 
        flash(form.username.errors[0], 'login_error')                                                
    if len(form.password.errors)>0:     
       flash(form.password.errors[0], 'login_error')    
    return redirect(request.values.get('next'))

def get_url_serializer(secret_key=None):
    if secret_key is None:
        secret_key = app.secret_key
    return URLSafeSerializer(secret_key)
    
@app.route('/register', methods=['GET', 'POST'])
def register():   
    if g.user is not None and g.user.is_authenticated():
        return redirect(url_for("index")) 
    if request.method == 'POST':        
        form = RegistrationForm(request.form)                
        if form.validate():            
            user = User()
            form.populate_obj(user)            
            user_exist = User.query.filter_by(username=form.username.data).first()
            mail_exist = User.query.filter_by(mail=form.mail.data).first()
            if user_exist:
                form.username.errors.append('Username already taken')
            if mail_exist:
                form.mail.errors.append('Email already in use')
            if user_exist or mail_exist:
                return flask.render_template('register.html', 
                                                register_form = form, 
                                                reset_form = ResetForm(), 
                                                login_form = LoginForm(), 
                                                active={}, 
                                                google_login_urls=google_login_urls(),
                                                notifications = get_notifications())
            else:                     
                now = datetime.datetime.now(timezone('UTC'))   
                user.date_joined = now              
                user.last_session = user.date_joined
                user.active = False 
                user.plan = "freemium"     
                user.password = user.hash_password(user.password)
                #validation_link
                s = get_url_serializer()
                validation_link = url_for('activate_user', payload=s.dumps(user.id), _external=True)
                msg = Message(
                  sender="no-reply@yarilab.com",
                  recipients=["%s"%user.mail],
                  html = mail_msg_body%(user.username, user.mail, user.mail, validation_link, validation_link) 
                )
                #mail.send(msg) 
                #TODO send mail                 
                db.session.add(user)                                                
                db.session.commit() 
                # create usage
                usage = Usage() 
                usage.user_id = user.id              
                usage.plan = "freemium"                
                usage.date = now.date()
                usage.queries = 0
                db.session.add(usage)                
                # create app
                mobile = Mobile()
                mobile.user_id = user.id
                os.makedirs(app.config['DATA_FOLDER']+'/%s/mobile'%user.id)
                os.makedirs(app.config['DATA_FOLDER']+'/%s/ivf'%user.id)                                 
                db.session.add(mobile)
                db.session.commit() #commit for usage, app
                #
                login_user(user)
                return flask.render_template('signup-success.html', 
                                            user = user, 
                                            active={}, 
                                            reset_form = ResetForm(), 
                                            login_form = LoginForm(),
                                            notifications = get_notifications())
        else:
            return flask.render_template('register.html', 
                                            register_form = form, 
                                            reset_form = ResetForm(), 
                                            login_form = LoginForm(), 
                                            active={}, 
                                            google_login_urls=google_login_urls(),
                                            notifications = get_notifications())      
    return flask.render_template('register.html', 
                                    register_form = RegistrationForm(), 
                                    reset_form = ResetForm(), 
                                    login_form = LoginForm(), 
                                    active={}, 
                                    google_login_urls=google_login_urls(),
                                    notifications = get_notifications())

@app.route('/user/activate/<payload>')
def activate_user(payload):
    s = get_url_serializer()
    try:
        user_id = s.loads(payload)
    except BadSignature:
        abort(404)
    user = User.query.get_or_404(user_id)
    user.active = True
    db.session.commit() 
    flash('User activated')
    return redirect(url_for('index'))

@app.route('/user/reset/<payload>')
def reset_password(payload):
    s = get_url_serializer()
    try:
        user_id = s.loads(payload)
    except BadSignature:
        abort(404)
    user = User.query.get_or_404(user_id)
    user.active = True
    db.session.commit() 
    flash('User activated')
    return redirect(url_for('index'))

        
@app.route("/logout")
@login_required
def logout():
    logout_user()
    session.clear() #TODO I added it laster after google login see if its ok ?
    return redirect(url_for("index"))

@app.route('/forgot_password', methods=['POST'])
def forgot_password():
    if g.user is not None and g.user.is_authenticated():
        return redirect(request.values.get('next'))
    reset_form = ResetForm(request.form)        
    if reset_form.validate_on_submit():
        user = User.query.filter_by(mail=reset_form.mail.data).first()                        
        if user:
            #TODO send mail
            flash('Reset link sent. Please check your inbox', 'reset_form_log')    
            return redirect(request.values.get('next'))
        else:
            flash('Please put a valid mail adress.', 'reset_form_error')                
            return redirect(request.values.get('next'))
    flash('Please put a valid mail adress', 'reset_form_error')    
    return redirect(request.values.get('next'))    
                                 
@app.route('/user/edit', methods=['GET', 'POST'])
@login_required
def user():
    """ user edit info """
    now = datetime.datetime.now()
    usage = Usage.query.filter(and_(Usage.user_id==current_user.id, func.YEAR(Usage.date)==now.year, func.MONTH(Usage.date)==now.month)).first()      
    if request.method == 'POST':
        edit_user_info_form = EditUserInfoForm(request.form)    
        if edit_user_info_form.validate_on_submit():
            current_user.username = edit_user_info_form.username.data
            current_user.mail = edit_user_info_form.mail.data            
            current_user.country = edit_user_info_form.country.data                        
            db.session.commit()                           
            return flask.render_template('user.html',
                     active={},
                     login_form = None,
                     reset_form = None,
                     edit_user_info_form = edit_user_info_form,
                     change_notif = 'success',
                     notifications = get_notifications())
        else:
            return flask.render_template('user.html',
                     active={},
                     login_form = None,
                     reset_form = None,
                     edit_user_info_form = edit_user_info_form,
                     paln = usage.plan,
                     change_notif = 'error',
                     notifications = get_notifications())                                                                
    return flask.render_template('user.html',
             active={},
             login_form = None, 
             reset_form = None,
             edit_user_info_form = EditUserInfoForm(),
             plan = usage.plan,
             change_notif = '',
             notifications = get_notifications())

@app.route('/user/password/change', methods=['GET', 'POST'])
@login_required
def change_password():
    """ user change password """
    if request.method == 'POST':
        change_pass_form = ChangePassForm(request.form)    
        if change_pass_form.validate_on_submit():    
            if current_user.verify_password(change_pass_form.old_password.data):        
                current_user.password = current_user.hash_password(change_pass_form.new_password.data)      
                db.session.commit()                                                   
                return flask.render_template('change.html',
                         active={},
                         login_form = None,
                         reset_form = None,
                         change_pass_form = change_pass_form,
                         notifications = get_notifications())                
            else:
                change_pass_form.old_password.errors.append('Wrong old password') 
                return flask.render_template('change.html',
                         active={},
                         login_form = None,
                         reset_form = None,
                         change_pass_form = change_pass_form,
                         notifications = get_notifications())                                     
        else:            
            return flask.render_template('change.html',
                     active={},
                     login_form = None,
                     reset_form = None,
                     change_pass_form = change_pass_form,
                     notifications = get_notifications())                                                              
    return flask.render_template('change.html',
                     active={},
                     login_form = None,
                     reset_form = None,
                     change_pass_form = ChangePassForm(),
                     notifications = get_notifications())

@app.route('/manage/projects/<query>', methods=['GET'])
def manage_database(query):
    if current_user.is_authenticated():                
        if query == "all":
            pps = Project.query.filter_by(user_id=current_user.id).all()                
        else:
            pps = Project.query.filter_by(user_id=current_user.id).whoosh_search(query, app.config.get('MAX_SEARCH_RESULTS')).all()     
        
        # usage details
        now = datetime.datetime.now()                
        usage = Usage.query.filter(and_(Usage.user_id==current_user.id, func.YEAR(Usage.date)==now.year, func.MONTH(Usage.date)==now.month)).first()
        api_usage = {"current_month":now.strftime("%B"), "visuals":{"now":current_user.visuals.count(), "limit":plans_limit[usage.plan]["visuals"]}, "queries":{"now":usage.queries , "limit":plans_limit[usage.plan]["queries"]}}              
                         
        # projects details 
        projects = {}
        for i,pp in enumerate(pps):
            #TODO update creation date to country timezone
            nb_visuals = pp.visuals.count() 
            nb_queries = Visual.query.filter_by(project_id=pp.id).join(Query).filter(Query.visual_id == Visual.id).count()
            projects[i] = {"project":pp, "nb_visuals":nb_visuals, "nb_queries":nb_queries}                            
        return flask.render_template('manage.html',
                 active={'manage':'class=active'},
                 login_form = None,
                 reset_form = None,
                 projects = projects,
                 project_search_form = SearchForm(),
                 api_usage = api_usage,
                 notifications = get_notifications())                                          
    """else:
        return flask.render_template('update.html',
                 active={'manage':'class=active'},
                 login_form = LoginForm(),
                 reset_form = ResetForm(),
                 project_search_form = SearchForm(),
                 google_login_urls=google_login_urls(),
                 notifications = get_notifications()) """
                 
@app.route('/manage/project/add', methods=['GET', 'POST'])
@login_required
def add_new_project():
    """ new project """
    if request.method == 'POST':
        new_project_form = NewProjectForm(request.form)    
        if new_project_form.validate_on_submit():
            project = Project()
            new_project_form.populate_obj(project)            
            project_name_exist = Project.query.filter_by(name=new_project_form.name.data, user_id=current_user.id).first() #TODO check double condition
            if project_name_exist:
                new_project_form.name.errors.append('project name already in use')
                return flask.render_template('new_project.html',
                    active={'manage':'class=active'},
                    login_form = None,
                    reset_form = None,
                    new_project_form = new_project_form,
                    notifications = get_notifications())
            project.creation_date = datetime.datetime.now(timezone('UTC'))               
            project.user = current_user
            db.session.add(project)
            db.session.commit()      
            project_dir = (app.config['DATA_FOLDER']+'/%s/%s')%(current_user.id, project.id)
            os.makedirs(project_dir+'/images')                 
            return redirect("/manage/projects/all")        
        else:
            return flask.render_template('new_project.html',
                    active={'manage':'class=active'},
                    login_form = None,
                    reset_form = None,                    
                    new_project_form = new_project_form,
                    notifications = get_notifications())        
    return flask.render_template('new_project.html',
                    active={'manage':'class=active'},
                    login_form = None,
                    reset_form = None,                    
                    new_project_form = NewProjectForm(),
                    notifications = get_notifications())

@app.route('/manage/project/edit/<project_id>', methods=['GET', 'POST'])
@login_required
def edit_project(project_id):
    """ edit project """
    if request.method == 'POST':
        edit_project_form = EditProjectForm(request.form)    
        if edit_project_form.validate_on_submit():
            project = Project.query.get(project_id)                    
            edit_project_form.populate_obj(project)            
            #update project
            db.session.commit() 
            return redirect("/manage/projects/all")        
        else:
            return flask.render_template('edit_project.html',
                    active={'manage':'class=active'},
                    login_form = None,
                    reset_form = None,                    
                    edit_project_form = edit_project_form,
                    project = project,
                    notifications = get_notifications())        
    #GET                
    project = Project.query.get(project_id)                
    return flask.render_template('edit_project.html',
                    active={'manage':'class=active'},
                    login_form = None,
                    reset_form = None,                    
                    edit_project_form = EditProjectForm(),
                    project = project,
                    notifications = get_notifications())
                                                                            
@app.route('/manage/project/delete/<project_id>', methods=['GET'])
@login_required
def delete_project(project_id):
    """ delete project: delete ivf index/images, visuals db, queries db, project db """
    project = Project.query.get(project_id)
    if project and project.user_id == current_user.id:        
        # delete physical data
        # TODO see problems of deleting when reading ?
        project_dir = (app.config['DATA_FOLDER']+'/%s/%s')%(current_user.id, project.id)
        if os.path.exists(project_dir):
            shutil.rmtree(project_dir)        
        for visual in project.visuals.all(): 
            for query in visual.queries.all():
                db.session.delete(query)                
            db.session.delete(visual)    
        db.session.delete(project)    
        db.session.commit()        
    return redirect("/manage/projects/all")                               
    
@app.route('/manage/visuals/edit/<project_id>/<query>', methods=['GET', 'POST'])
@login_required
def edit_project_visuals(project_id, query):
    """ edit project's visuals """
    project = Project.query.get(project_id)
    if project and project.user_id == current_user.id:        
        project_name = project.name
        #total visuals
        total = Visual.query.filter_by(project_id = project_id).count()    
        # get page items
        page = int(request.args.get('page', 1))
        per_page = request.args.get('per_page')
        if not per_page:    per_page = app.config.get('PER_PAGE', 10)
        else:   per_page = int(per_page)
        offset = (page - 1) * per_page    
        # get visuals subset
        if query=="all":
            visuals = Visual.query.filter_by(project_id = project_id)\
                                        .order_by(Visual.creation_date)\
                                        .offset(offset)\
                                        .limit(per_page)\
                                        .all()                                            
        else:
            #TODO manage visuals query search
            #visuals = Visual.query.filter_by(user_id=current_user.id).whoosh_search(query, app.config.get('MAX_SEARCH_RESULTS')).all()  
            visuals = Visual.query.filter_by(project_id = project_id)\
                                        .order_by(Visual.creation_date)\
                                        .offset(offset)\
                                        .limit(per_page)\
                                        .all()                                            
                                                                                    
        pagination = Pagination(page=page, 
                                per_page=per_page, 
                                total=total, 
                                record_name='visuals', 
                                css_framework=app.config.get('CSS_FRAMEWORK', 'bootstrap3'), 
                                link_size=app.config.get('LINK_SIZE', 'sm'), 
                                show_single_page=app.config.get('SHOW_SINGLE_PAGE', False)
                                )        
        return flask.render_template('edit_visuals.html',
                    active={'manage':'class=active'},
                    login_form = None,
                    reset_form = None,
                    visuals = visuals,
                    page = page,
                    per_page = per_page,
                    pagination = pagination,                
                    project_id = project_id,
                    project_name = project_name,
                    notifications = get_notifications())    
    else:
        return redirect(url_for('index'))
        
@app.route('/manage/visual/edit/<visual_id>', methods=['GET', 'POST'])
@login_required
def edit_visual(visual_id):
    """ edit single visual """          
    visual = Visual.query.get(visual_id)
    project = Project.query.get(visual.project_id)
    project_name = project.name
    if request.method == 'POST' and project and project.user_id == current_user.id:
        edit_visual_form = EditVisualForm(request.form)                        
        if edit_visual_form.validate():
            edit_visual_form.populate_obj(visual)            
            visual.creation_date = datetime.datetime.now(timezone('UTC'))
            db.session.commit()
            return redirect('/manage/visuals/edit/%s'%visual.project_id)         
        else:
            return flask.render_template('edit_visual.html',
                        active={'manage':'class=active'},
                        login_form = None,
                        reset_form = None,
                        visual = visual,
                        project_name = project_name,
                        edit_visual_form = edit_visual_form,
                        notifications = get_notifications())                                   
    return flask.render_template('edit_visual.html',
                active={'manage':'class=active'},
                login_form = None,
                reset_form = None,
                visual = visual,
                project_name = project_name,
                edit_visual_form = EditVisualForm(),
                notifications = get_notifications())        

"""@app.route('/manage/<project_id>/visual/add/', methods=['GET', 'POST'])
@login_required
def add_visual(project_id):
    #add single visual 
    print "req:", request.method        
    project = Project.query.get(project_id)
    project_name = project.name    
    if request.method == 'POST' and project and project.user_id == current_user.id:        
        add_visual_form = EditVisualForm(request.form)                        
        if add_visual_form.validate():
            visual = Visual()
            add_visual_form.populate_obj(visual)            
            visual.creation_date = datetime.datetime.now(timezone('UTC'))
            db.session.commit()
            return redirect('/manage/visuals/edit/%s'%visual.project_id)         
        else:
            return flask.render_template('add_visual.html',
                        active={'manage':'class=active'},
                        login_form = None,
                        reset_form = None,
                        project_id = project_id,
                        project_name = project_name,
                        add_visual_form = add_visual_form,
                        upload_visual_error = ""
                        )                                   
    return flask.render_template('add_visual.html',
                active={'manage':'class=active'},
                login_form = None,
                reset_form = None,
                project_id = project_id,
                project_name = project_name,
                add_visual_form = EditVisualForm(),
                upload_visual_error = ""                
                )        
"""
# ------------------------------------------------ UPLOAD STUFF ---------------------------------------------- #
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS           

def get_file_size(file):
    file.seek(0, 2)  # Seek to the end of the file
    size = file.tell()  # Get the position of EOF
    file.seek(0)  # Reset the file position to the beginning
    return size

@app.route('/manage/<project_id>/visuals/add/', methods=['GET', 'POST', 'DELETE'])
@login_required
def add_visuals(project_id):
    """ add visuals """
    # TODO see: https://pythonhosted.org/Flask-Uploads/ is better for uploads handling
    project = Project.query.get(project_id)
    project_name = project.name      
    if request.method == 'POST' and project and project.user_id == current_user.id:                                
        results = []        
        for key, file in request.files.iteritems():
            if file:
                result = {}
                name = secure_filename(file.filename)                
                size = get_file_size(file)
                result['name'] = name 
                result['size'] = size            
                if allowed_file(file.filename):
                    usage = Usage.query.filter(and_(Usage.user_id==current_user.id, func.YEAR(Usage.date)==now.year, func.MONTH(Usage.date)==now.month)).first()
                    if current_user.visuals.count() >= plans_limit[usage.plan]['visuals'] :                    
                        #TODO manage exceeded visuals 
                        result['error'] = "Visual not uploaded, your plan's limit is reached. Please upgrade your plan." 
                        results.append(result)                
                        return jsonify(files=results)  
                    #file informations
                    hashed_name = hashlib.md5(datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S-%f")).hexdigest()+'_'+hashlib.md5('.'.join(name.split('.')[0:-1])).hexdigest()+'.'+name.split('.')[-1]
                    url = os.path.join((app.config['DATA_FOLDER']+'/%s/%s/images/')%(current_user.id, project_id), hashed_name)                           
                    meta = request.form['meta']
                    visual = Visual()
                    visual.user_id = current_user.id
                    visual.name = name
                    visual.visual_url = url
                    visual.visual_meta = meta
                    visual.project_id = project_id
                    visual.creation_date = datetime.datetime.now(timezone('UTC'))                                        
                    file.save(url)                                                            
                    db.session.add(visual)                    
                    db.session.commit()                      
                    #Image.open(image).thumbnail(size).save("thumbnail_%s" % image)          
                    result['url'] = url
                    result['delete_type'] = 'DELETE'
                    result['delete_url'] = '/visual/delete/%s'%visual.id
                    result['thumbnailUrl'] = url
                    result['meta'] = meta                     
                    result['id'] = visual.id
                else:
                    result['error'] = 'Filetype not allowed'                
                results.append(result)                
        return jsonify(files=results)                                            
    return flask.render_template('add_visuals.html',
                active={'manage':'class=active'},
                login_form = None,
                reset_form = None,
                project_id = project_id,
                project_name = project_name ,
                notifications = get_notifications())  

@app.route('/visual/delete/<int:visual_id>', methods=['DELETE'])
@login_required
def delete_visual(visual_id):
    """ delete visual image """
    visual = Visual.query.get(visual_id)        
    project = Project.query.get(visual.project_id)     
    if (request.method == 'DELETE') and (project.user_id == current_user.id):                                        
        # TODO remove visual index from ivf
        os.remove(visual.visual_url)
        db.session.delete(visual)
        db.session.commit()   
        return jsonify({"status":"success"})
    return jsonify({"status":"denied"})  
                   
# ------------------------------------------------ ./UPLOAD STUFF ---------------------------------------------- #    
    
"""@app.route('/data/<path:filename>')
def custom_static(filename):
    return flask.send_from_directory("."+app.config['DATA_FOLDER']+"/", filename)
"""            

@app.route('/dashboard', methods=["GET"])
@login_required
def dashboard():
    projects = Project.query.filter_by(user_id = current_user.id).all()    
    return flask.render_template('dashboard.html',
             active={'manage': 'class=active'},
             login_form = None,
             reset_form = None,
             projects = projects,
             notifications = get_notifications())

@app.route('/dashboard/<int:project_id>', methods=["POST"])
@login_required
def dashboard_series(project_id):
    if request.method == 'POST':        
        action = request.form.get('action')
        #years = map(lambda x:x.quering_date.year, Query.query.join(Visual).filter(Visual.project_id == project_id).group_by(func.year(Query.quering_date)).all())
        series = {} 
        unique_users = {}
        if action == "get_queries":                            
            queries = Query.query.join(Visual).filter(Visual.project_id == project_id).all()                        
            if len(queries) > 0:                
                center_markers = {}
                for query in queries:
                    year = query.quering_date.year
                    month = query.quering_date.month
                    day = query.quering_date.day
                    query_user = query.phone_id    
                    latitude = query.latitude
                    longitude = query.longitude                    
                    city = query.city 
                    pos = '(%f, %f)'%(latitude, longitude)                                                     
                    if year in series.keys():
                        if query_user in series[year]['users_markers'].keys():
                            if pos in series[year]['users_markers'][query_user].keys():
                                series[year]['users_markers'][query_user][pos] += 1 #freq incrementation
                            else:
                                series[year]['users_markers'][query_user][pos] = 1  
                        else:     
                            series[year]['users_markers'][query_user] = {pos:1}
                        if query_user in unique_users[year].keys():    
                            unique_users[year][query_user] += 1                  
                        else:
                            unique_users[year][query_user] = 1                            
                        """if city in series[year]['cities']:
                            series[year]['cities'][city] += 1
                        else:         
                            series[year]['cities'][city] = 1"""
                        series[year]['details'][month]['details'][day] += 1
                        series[year]['details'][month]['overall'] += 1
                        series[year]['overall'] += 1                    
                    else: 
                        center_markers[year] = [0.0, 0.0]  
                        unique_users[year] = {}
                        series[year] = { 'overall':0, 'details':{1:{}, 2:{}, 3:{}, 4:{}, 5:{}, 6:{}, 7:{}, 8:{}, 9:{}, 10:{}, 11:{}, 12:{}}, 'users_markers':{}, 'center_markers':None, 'cities':{} }
                        for i in range(1, 13):
                            nb_days = monthrange(year, i)[1]
                            series[year]['details'][i]['overall'] = 0
                            series[year]['details'][i]['details'] = {}
                            for j in range(1, nb_days+1):
                                series[year]['details'][i]['details'][j] = 0                                
                        series[year]['details'][month]['details'][day] = 1               
                        series[year]['details'][month]['overall'] = 1                                   
                        series[year]['overall'] = 1                             
                        unique_users[year][query_user] = 1  
                        series[year]['users_markers'][query_user] = {pos:1}                
                        #series[year]['cities'][city] = 1
                    center_markers[year][0] += latitude
                    center_markers[year][1] += longitude    
                # users handler
                for year in series.keys():        
                    series[year]['unique_users'] = len(unique_users[year].keys())
                    #Remark: returning users are those who queried more than once, or more than 1 day ? for instance its more than 1 query
                    series[year]['returning_users'] = 0
                    for key, value in unique_users[year].iteritems():
                        if value > 1:
                            series[year]['returning_users'] += 1.0
                    #series[year]['returning_users'] = (series[year]['returning_users']/len(unique_users[year].keys()))*100                   
                    # for test:          
                    series[year]['returning_users'] = 55                        
                    # center markers
                    center_markers[year][0] /= series[year]['overall']
                    center_markers[year][1] /= series[year]['overall']
                    series[year]['center_markers'] = center_markers[year] 
                    # get cities + frequencies + sorted
                    cities = db.session.query(Query.city, func.count('*').label('qc'))\
                                .join(Visual)\
                                .filter(and_(Visual.project_id == project_id, func.Year(Query.quering_date) == year))\
                                .group_by(Query.city)\
                                .order_by(desc('qc'))\
                                .all()
                    ss = sum(map(lambda x:x[1], cities))            
                    cities = map(lambda x:x+('%.1f'%(x[1]*100.0/ss) ,), cities)
                    series[year]['cities'] = cities

            else:
                #zero queries handler
                year = datetime.datetime.now().year
                series[year] = {'overall':0, 'details':{1:{}, 2:{}, 3:{}, 4:{}, 5:{}, 6:{}, 7:{}, 8:{}, 9:{}, 10:{}, 11:{}, 12:{}}, 'returning_users':0, 'unique_users':0}        
                for i in range(1, 13):
                    series[year]['details'][i]['overall'] = 0
                    series[year]['details'][i]['details'] = {}                    
                    nb_days = monthrange(year, i)[1]
                    for j in range(1, nb_days+1):
                        series[year]['details'][i]['details'][j] = 0                                    
            return jsonify({"series": series})
        elif action == "rank_visuals":                
            series = db.session.query(Visual.visual_url, func.count('*').label('qc'))\
                        .join(Query).filter(Visual.project_id == project_id)\
                        .group_by(Query.visual_id)\
                        .order_by(desc('qc')).limit(10).all()
            return jsonify({"series": series})        
        else:
            return jsonify({"status":"denied"})              
    return jsonify({"status":"denied"})          
                                                                                               
@login_manager.request_loader
def load_user_from_request(request):
    # first, try to login using the api_key url arg
    api_key = request.args.get('api_key')
    if api_key:
        user = User.query.filter_by(api_key=api_key).first()
        if user:
            return user
    # next, try to login using Basic Auth
    api_key = request.headers.get('Authorization')
    if api_key:
        api_key = api_key.replace('Basic ', '', 1)
        try:
            api_key = base64.b64decode(api_key)
        except TypeError:
            pass
        user = User.query.filter_by(api_key=api_key).first()
        if user:
            return user
    # finally, return None if both methods did not login the user
    return None
        
@app.route('/api/users/<int:id>')
def get_user(id):
    user = User.query.get(id)
    if not user:
        abort(400)
    return jsonify({'username': user.username})
    
@app.route('/api/token')
@login_required
def get_auth_token():
    token = g.user.generate_auth_token(600) #in seconds
    return jsonify({'token': token.decode('ascii'), 'duration': 600})


@app.route('/api/resource')
@login_required
def get_resource():
    return jsonify({'data': 'Hello, %s!' % g.user.username})
    
@app.route('/api/add', methods=['POST'])                                     
@login_required
def add():
    try:
        # We will save the file to disk for possible data collection.
        imagefile = request.files['imagefile']
        filename = os.path.join(FLAGS.upload_folder,
                                str(datetime.datetime.now(timezone('UTC'))).replace(' ', '_') + \
                                secure_filename(imagefile.filename))
        imagefile.save(filename)
        logging.info('Saving uploaded to %s.', filename)        
    except Exception as err:
        error = 'Add upload image error: %s', err
        logging.info(error)
        return {"success":"nope", "error":error} 
    #app.model.add(filename)
    # TODO keep only thimbnail or remove filename ?
    return jsonify({"success":"added", "error":"nope"}) 


def check_usage():
    pass

def update_usage():
    # update query table
    # update usage table
    # update notification table if excess
    pass
      
@app.route('/cart/<requested_plan>')
@login_required
def cart(requested_plan):
    if requested_plan in ['basic', 'premium']:
        now = datetime.datetime.now()
        usage = Usage.query.filter(and_(Usage.user_id==current_user.id, func.YEAR(Usage.date)==now.year, func.MONTH(Usage.date)==now.month)).first()
        actual_plan = usage.plan
        return flask.render_template('cart.html',
                                     active={},
                                     login_form = None,
                                     reset_form = None,
                                     actual_plan = actual_plan,
                                     requested_plan = requested_plan,
                                     notifications = get_notifications())
        return flask.render_template('index.html',
                                     active={},
                                     login_form = None,
                                     reset_form = None,
                                     notifications = get_notifications())

@app.route('/query')
@login_required
def query():  
    return flask.render_template('query.html',
                                 active={'manage':'class=active'},
                                 login_form = None,
                                 reset_form = None,
                                 mobile = current_user.mobile.first(),
                                 notifications = get_notifications())

@app.route('/query/send', methods=['POST'])
@login_required
def send_link():
    # compile app
    mobile = current_user.mobile
    # create tmp project file + launch compile cmd
    
    # send link by mail
    s = get_url_serializer()
    validation_link = url_for('activate_user', payload=s.dumps(user.id), _external=True)
    msg = Message(
      sender="no-reply@yarilab.com",
      recipients=["%s"%current_user.mail],
      html = mail_msg_body%(current_user.username, current_user.mail, current_user.mail, validation_link, validation_link) 
    )
    #mail.send(msg) 
    send_state = "success"
    return flask.render_template('query.html',
                                 active={'manage':'class=active'},
                                 login_form = None,
                                 reset_form = None,
                                 mobile = current_user.mobile.first(),
                                 send_state = send_state,
                                 notifications = get_notifications())

@app.route('/query/edit', methods=['GET', 'POST'])
@login_required
def edit_query():
    if request.method == 'POST':                
        edit_mobile_form = EditMobileForm(request.form)                        
        if edit_mobile_form.validate():
            mobile = current_user.mobile.first()            
            logo_url = mobile.logo_url            
            icon_url = mobile.icon_url
            edit_mobile_form.populate_obj(mobile)
            # form has formal log/icon url (fields just for error rendering)
            mobile.logo_url = logo_url
            mobile.icon_url = icon_url            
            key = request.files.keys()[0]
            file = request.files[key]
            if len(file.filename)>0:                
                if allowed_file(file.filename):
                    name = secure_filename(file.filename)  
                    hashed_name = hashlib.md5(datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S-%f")).hexdigest()+'_'+hashlib.md5('.'.join(name.split('.')[0:-1])).hexdigest()+'.'+name.split('.')[-1]    
                    url = os.path.join(app.config['DATA_FOLDER']+'/%s/mobile/'%current_user.id, hashed_name)                           
                    file.save(url) 
                    if key=="logo": 
                        if mobile.logo_url != "":
                            os.remove(mobile.logo_url)
                        mobile.logo_url = url
                    else:                           
                        if mobile.icon_url != "":                            
                            os.remove(mobile.icon_url)                                                    
                        mobile.icon_url = url
                else:
                    if key=="logo":
                        edit_mobile_form.logo_url.errors.append("not allowed file format")
                    else:
                        edit_mobile_form.icon_url.errors.append("not allowed file format")                            
                    return flask.render_template('edit_app.html',
                             active={'manage':'class=active'},
                             login_form = None,
                             reset_form = None,                                 
                             edit_mobile_form = edit_mobile_form,
                             notifications = get_notifications())                                                         
            db.session.commit()  #we save other parameters anyway            
            return redirect('/query')
        else:   
            return flask.render_template('edit_app.html',
                     active={'manage':'class=active'},
                     login_form = None,
                     reset_form = None,                                 
                     edit_mobile_form = edit_mobile_form,
                     notifications = get_notifications())            
    edit_mobile_form = EditMobileForm(obj=current_user.mobile.first())                      
    return flask.render_template('edit_app.html',
                     active={'manage':'class=active'},
                     login_form = None,
                     reset_form = None,                                 
                     edit_mobile_form = edit_mobile_form,
                     notifications = get_notifications())

@app.route('/user/usage')
@login_required
def usage():
    now = datetime.datetime.now()    
    usage = Usage.query.filter(and_(Usage.user_id==current_user.id, func.Month(Usage.date)==now.month, func.Year(Usage.date)==now.year )).first()
    api_usage = {"plan":usage.plan, "plan_expiring_date":usage.plan_expiring_date, "current_month":now.strftime("%B"), "visuals":{"now":current_user.visuals.count(), "limit":plans_limit[usage.plan]["visuals"]}, "queries":{"now":usage.queries , "limit":plans_limit[usage.plan]["queries"], "query_excess_price":plans_limit[usage.plan]["query_excess_price"]}}  
    bills = Billing.query.filter(Billing.user_id==current_user.id).order_by(desc(Billing.paiement_date)).all()
    return flask.render_template('usage.html',
             active={},
             login_form = None,
             reset_form = None,
             google_login_urls = None,
             bills = bills,
             api_usage = api_usage,
             notifications = get_notifications())
                                                                                                                                                                                                                                                        
@app.route('/about')
def about():
    return flask.render_template('about.html',
                                 active={'about':'class=active'},
                                 login_form = LoginForm(),
                                 reset_form = ResetForm(),
                                 google_login_urls = google_login_urls(),
                                 notifications = get_notifications())

@app.route('/notifications/defresh', methods=['POST'])
@login_required
def defresh():
    notifs = Notification\
            .query.filter(and_(Notification.user_id == current_user.id, Notification.fresh == True))\
            .update(dict(fresh = False))
    db.session.commit()
    return jsonify(response={"state":"ok"})

@app.route('/notifications/all', methods=['POST'])
@login_required
def all_notifications():    
    notifications = get_notifications(limit="all")
    result = []
    for notification in notifications:
        result.append({"post":notification.post, "creation_date":notification.creation_date, "notif_type":notification.notif_type})
    return jsonify(notifications = result)  
          
    
@app.errorhandler(404)
def not_found(error):
    return flask.render_template('index.html', login_form = LoginForm(), reset_form = ResetForm(), active={}, google_login_urls = google_login_urls()), 404


@app.route('/')
def index():        
    return flask.render_template('index.html',
                                 active={},
                                 login_form = LoginForm(),
                                 reset_form = ResetForm(),
                                 google_login_urls = google_login_urls(),
                                 notifications = get_notifications())

@app.route('/demo')
def demo():  
    return flask.render_template('demo.html',
                                 active={'demo':'class=active'},
                                 login_form = LoginForm(),
                                 reset_form = ResetForm(),
                                 google_login_urls = google_login_urls(),
                                 notifications = get_notifications())

@app.route('/offer')
def offer():  
    return flask.render_template('offer.html',
                                 active={'offer':'class=active'},
                                 login_form = LoginForm(),
                                 reset_form = ResetForm(),
                                 google_login_urls = google_login_urls(),
                                 notifications = get_notifications())
                                 
@app.route('/pricing')
def pricing():
    return flask.render_template('pricing.html',
                                 active={'pricing':'class=active'},
                                 login_form = LoginForm(),
                                 reset_form = ResetForm(),
                                 google_login_urls = google_login_urls(),
                                 notifications = get_notifications())


app.route('/match2', methods=['POST'])
def match2():
    return jsonify({"denied":"Access denied to DB"})	


@app.route('/', subdomain="api", methods=['POST'])                                     
@login_required
def match():
    # usage checking and updating
    API_state = check_usage()
    if API_state:                        
        # classify image using the image name
        try:
            # We will save the file to disk for possible data collection.
            imagefile = request.files['imagefile']
            #pdb.set_trace()
            filename = os.path.join(FLAGS.upload_folder,
                                    str(datetime.datetime.now(timezone('UTC'))).replace(' ', '_') + \
                                    secure_filename(imagefile.filename))
            imagefile.save(filename)
            logging.info('  Saving uploaded to %s', filename)        
        except Exception as err:
            error = 'Match upload image got error: %s', err
            logging.info(error)
            return {"success":"nope", "error":error} 
        #matches = app.model.match(filename)
        # TODO keep only thimbnail or remove filename ?
        return jsonify({"success":"matched", "error":"nope"})
    return jsonify({"denied":"Access denied to DB"})          

                                     
