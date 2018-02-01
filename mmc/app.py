# -*- coding: utf_8 -*-

import os

from flask import Flask
from flask import render_template
from flask import url_for

from mmc.modules.api.views import api
from mmc.modules.page.views import page
from mmc.extensions import debug_toolbar
from mmc.utils import ReverseProxied


def create_app(settings_override=None):
    """
    Create a Flask application using the app factory pattern.

    :param settings_override: Override settings
    :return: Flask app
    """
    app = Flask(__name__, instance_relative_config=True)

    app.config.from_object('config.settings')

    if app.config.from_envvar('MMCVIEWER_SETTINGS', silent=True):
        env_settings = os.environ['MMCVIEWER_SETTINGS']
        app.logger.info('Using MMCVIEWER_SETTINGS: {}'.format(env_settings))

    if settings_override:
        app.logger.info('Overriding settings with parameters')
        app.config.update(settings_override)

    app.logger.setLevel(app.config['LOG_LEVEL'])

    middleware(app)

    app.register_blueprint(api)
    app.register_blueprint(page)

    extensions(app)
    error_templates(app)

    return app


def extensions(app):
    """
    Register 0 or more extensions (mutates the app passed in).

    :param app: Flask application instance
    :return: None
    """
    debug_toolbar.init_app(app)

    return None


def middleware(app):
    """
    Register 0 or more middleware (mutates app that is passed in).

    :param app: Flask application instance
    :return: None
    """
    app.wsgi_app = ReverseProxied(app.wsgi_app)

    return None


def error_templates(app):
    """
    Register 0 or more custom error pages (mutates the app passed in).

    :param app: Flask application instance
    :return: None
    """

    def render_status(status):
        """
         Render a custom template for a specific status.
           Source: http://stackoverflow.com/a/30108946

         :param status: Status as a written name
         :type status: str
         :return: None
         """
        # Get the status code from the status, default to a 500 so that we
        # catch all types of errors and treat them as a 500.
        code = getattr(status, 'code', 500)
        redirect_url = url_for('page.index')

        # for specific error pages, you could create an errors/404.html or an
        # errors/500.html and than do something like the following
        #
        # return render_template('errors/{0}.html'.format(code),
        #                        redirect_url=redirect_url), code

        return render_template('errors/error.html',
                               error_code=code,
                               redirect_url=redirect_url), code

    for error in [404, 500]:
        app.errorhandler(error)(render_status)

    return None

